import openmm as mm
import openmm.unit as unit
from openmm.app import PDBFile, PDBxFile
from pdbfixer.pdbfixer import PDBFixer, proteinResidues, dnaResidues, rnaResidues, _guessFileFormat
from flask import Flask, request, session, g, render_template, make_response, send_file, url_for
from werkzeug.utils import secure_filename
from werkzeug.serving import make_server
from multiprocessing import Process, Pipe
import datetime
from io import StringIO
import argparse
import os
import shutil
import signal
import sys
import tempfile
import threading
import time
import traceback
import webbrowser
import zipfile

app = Flask(__name__)
app.config.from_object(__name__)
app.config.update({'SECRET_KEY':'development key'})
app.jinja_env.globals['mm'] = mm

uploadedFiles = {}
fixer = None
scriptOutput = None
simulationProcess = None

def saveUploadedFiles():
    uploadedFiles.clear()
    for key in request.files:
        filelist = []
        for file in request.files.getlist(key):
            temp = tempfile.TemporaryFile()
            shutil.copyfileobj(file, temp)
            filelist.append((temp, secure_filename(file.filename)))
        uploadedFiles[key] = filelist

@app.route('/headerControls')
def headerControls():
    if 'startOver' in request.args:
        return showSelectFileType()
    if 'quit' in request.args:
        shutdownEvent.set()
        return "OpenMM Setup has stopped running.  You can close this window."

@app.route('/')
def showSelectFileType():
    return render_template('selectFileType.html')

@app.route('/selectFiles')
def selectFiles():
    session['fileType'] = request.args.get('type', '')
    return showConfigureFiles()

def showConfigureFiles():
    try:
        fileType = session['fileType']
        if fileType == 'pdb':
            return render_template('configurePdbFile.html')
        elif fileType == 'amber':
            return render_template('configureAmberFiles.html')
        elif fileType == 'charmm':
            return render_template('configureCharmmFiles.html')
        elif fileType == 'gromacs':
            return render_template('configureGromacsFiles.html')
        elif fileType == 'tinker':
            return render_template('configureTinkerFiles.html')
    except:
        app.logger.error('Error displaying configure files page', exc_info=True)
    # The file type is invalid, so send them back to the select file type page.
    return showSelectFileType()

@app.route('/configureFiles', methods=['POST'])
def configureFiles():
    fileType = session['fileType']
    if fileType == 'pdb':
        if 'file' not in request.files or request.files['file'].filename == '':
            # They didn't select a file.  Send them back.
            return showConfigureFiles()
        saveUploadedFiles()
        session['forcefield'] = request.form.get('forcefield', '')
        session['waterModel'] = request.form.get('waterModel', '')
        session['cleanup'] = request.form.get('cleanup', '')
        file, name = uploadedFiles['file'][0]
        file.seek(0, 0)
        session['pdbType'] = _guessFileFormat(file, name)
        if session['cleanup'] == 'yes':
            global fixer
            if session['pdbType'] == 'pdb':
                fixer = PDBFixer(pdbfile=file)
            else:
                fixer = PDBFixer(pdbxfile=StringIO(file.read().decode()))
            configureDefaultOptions()
            return showSelectChains()
    elif fileType == 'amber':
        if 'prmtopFile' not in request.files or request.files['prmtopFile'].filename == '' or 'inpcrdFile' not in request.files or request.files['inpcrdFile'].filename == '':
            # They didn't select a file.  Send them back.
            return showConfigureFiles()
        saveUploadedFiles()
    elif fileType == 'charmm':
        if 'psfFile' not in request.files or request.files['psfFile'].filename == '' or 'crdFile' not in request.files or request.files['crdFile'].filename == '':
            # They didn't select a file.  Send them back.
            return showConfigureFiles()
        saveUploadedFiles()
    elif fileType == 'gromacs':
        if 'topFile' not in request.files or request.files['topFile'].filename == '' or 'groFile' not in request.files or request.files['groFile'].filename == '':
            # They didn't select a file.  Send them back.
            return showConfigureFiles()
        saveUploadedFiles()
        session['gromacsIncludeDir'] = request.form.get('gromacsIncludeDir', '')
    elif fileType == 'tinker':
        if 'xyzFile' not in request.files or request.files['xyzFile'].filename == '':
            # They didn't select a file.  Send them back.
            return showConfigureFiles()
        session['waterModel'] = request.form.get('waterModel', '')
        saveUploadedFiles()
    configureDefaultOptions()
    return showSimulationOptions()

@app.route('/getCurrentStructure')
def getCurrentStructure():
    pdb = StringIO()
    PDBFile.writeFile(fixer.topology, fixer.positions, pdb)
    return pdb.getvalue()

def showSelectChains():
    chains = []
    hasHeterogen = False
    for chain in fixer.topology.chains():
        residues = list(r.name for r in chain.residues())
        if any(r in proteinResidues for r in residues):
            content = "Protein"
        elif any(r in rnaResidues for r in residues):
            content = "RNA"
        elif any(r in dnaResidues for r in residues):
            content = "DNA"
        else:
            content = ', '.join(set(residues))
            hasHeterogen = True
        chains.append((chain.id, len(residues), content))
    if len(chains) < 2 and not hasHeterogen:
        session['heterogens'] = 'all'
        return showAddResidues()
    return render_template('selectChains.html', chains=chains)

@app.route('/selectChains', methods=['POST'])
def selectChains():
    session['heterogens'] = request.form.get('heterogens', '')
    numChains = len(list(fixer.topology.chains()))
    request.form.getlist('include')
    deleteIndices = [i for i in range(numChains) if str(i) not in request.form.getlist('include')]
    fixer.removeChains(deleteIndices)
    return showAddResidues()

def showAddResidues():
    spans = []
    chains = list(fixer.topology.chains())
    fixer.findMissingResidues()
    if len(fixer.missingResidues) == 0:
        return showConvertResidues()
    for i, key in enumerate(sorted(fixer.missingResidues)):
        residues = fixer.missingResidues[key]
        chain = chains[key[0]]
        chainResidues = list(chain.residues())
        if key[1] < len(chainResidues):
            offset = int(chainResidues[key[1]].id)-len(residues)-1
        else:
            offset = int(chainResidues[-1].id)
        spans.append((chain.id, offset+1, offset+len(residues), ', '.join(residues)))
    return render_template('addResidues.html', spans=spans)

@app.route('/addResidues', methods=['POST'])
def addResidues():
    keys = [key for key in sorted(fixer.missingResidues)]
    for i, key in enumerate(keys):
        if str(i) not in request.form.getlist('add'):
            del fixer.missingResidues[key]
    return showConvertResidues()

def showConvertResidues():
    fixer.findNonstandardResidues()
    if len(fixer.nonstandardResidues) == 0:
        return showAddHeavyAtoms()
    residues = []
    nucleotides = ['DA', 'DC', 'DG', 'DT', 'A', 'C', 'G', 'T']
    for i in range(len(fixer.nonstandardResidues)):
        residue, replaceWith = fixer.nonstandardResidues[i]
        if replaceWith in proteinResidues:
            replacements = proteinResidues
        else:
            replacements = nucleotides
        residues.append((residue.chain.id, residue.name, residue.id, replacements, replaceWith))
    return render_template('convertResidues.html', residues=residues)

@app.route('/convertResidues', methods=['POST'])
def convertResidues():
    for i in range(len(fixer.nonstandardResidues)):
        if str(i) in request.form.getlist('convert'):
            fixer.nonstandardResidues[i] = (fixer.nonstandardResidues[i][0], request.form['residue'+str(i)])
    fixer.replaceNonstandardResidues()
    return showAddHeavyAtoms()

def showAddHeavyAtoms():
    if session['heterogens'] == 'none':
        fixer.removeHeterogens(False)
    elif session['heterogens'] == 'water':
        fixer.removeHeterogens(True)
    fixer.findMissingAtoms()
    allResidues = list(set(fixer.missingAtoms.keys()).union(fixer.missingTerminals.keys()))
    allResidues.sort(key=lambda x: x.index)
    if len(allResidues) == 0:
        return addHeavyAtoms()
    residues = []
    for residue in allResidues:
        atoms = []
        if residue in fixer.missingAtoms:
            atoms.extend(atom.name for atom in fixer.missingAtoms[residue])
        if residue in fixer.missingTerminals:
            atoms.extend(atom for atom in fixer.missingTerminals[residue])
        residues.append((residue.chain.id, residue.name, residue.id, ', '.join(atoms)))
    return render_template('addHeavyAtoms.html', residues=residues)

@app.route('/addHeavyAtoms', methods=['POST'])
def addHeavyAtoms():
    fixer.addMissingAtoms()
    return showAddHydrogens()

def showAddHydrogens():
    unitCell = fixer.topology.getUnitCellDimensions()
    if unitCell is not None:
        unitCell = unitCell.value_in_unit(unit.nanometer)
    boundingBox = tuple((max((pos[i] for pos in fixer.positions))-min((pos[i] for pos in fixer.positions))).value_in_unit(unit.nanometer) for i in range(3))
    return render_template('addHydrogens.html', unitCell=unitCell, boundingBox=boundingBox)

@app.route('/addHydrogens', methods=['POST'])
def addHydrogens():
    if 'addHydrogens' in request.form:
        pH = float(request.form.get('ph', '7'))
        fixer.addMissingHydrogens(pH)
    if 'addWater' in request.form:
        padding, boxSize, boxShape = None, None, None
        if request.form['boxType'] == 'geometry':
            padding = float(request.form['geomPadding']) * unit.nanometer
            boxShape = request.form['geometryDropdown']
        else:
            boxSize = (float(request.form['boxx']), float(request.form['boxy']), float(request.form['boxz']))*unit.nanometer
        ionicStrength = float(request.form['ionicstrength'])*unit.molar
        positiveIon = request.form['positiveion']+'+'
        negativeIon = request.form['negativeion']+'-'
        fixer.addSolvent(boxSize=boxSize, padding=padding, boxShape=boxShape, positiveIon=positiveIon, negativeIon=negativeIon, ionicStrength=ionicStrength)
    elif 'addMembrane' in request.form:
        lipidType = request.form['lipidType']
        padding = float(request.form['membranePadding'])*unit.nanometer
        ionicStrength = float(request.form['ionicstrength'])*unit.molar
        positiveIon = request.form['positiveion']+'+'
        negativeIon = request.form['negativeion']+'-'
        fixer.addMembrane(lipidType=lipidType, minimumPadding=padding, positiveIon=positiveIon, negativeIon=negativeIon, ionicStrength=ionicStrength)
    
    # Save the new PDB file.
    
    uploadedFiles['originalFile'] = uploadedFiles['file']
    pdb = StringIO()
    if session['pdbType'] == 'pdb':
        try:
            PDBFile.writeFile(fixer.topology, fixer.positions, pdb, True)
        except:
            # This can happen if the ids are too large to fit in the allowed space.
            pdb = StringIO()
            PDBFile.writeFile(fixer.topology, fixer.positions, pdb, False)
    else:
        PDBxFile.writeFile(fixer.topology, fixer.positions, pdb, True)
    temp = tempfile.TemporaryFile()
    temp.write(pdb.getvalue().encode('utf-8'))
    name = uploadedFiles['file'][0][1]
    dotIndex = name.rfind('.')
    if dotIndex == -1:
        prefix = name
        suffix = ''
    else:
        prefix = name[:dotIndex]
        suffix = name[dotIndex:]
    uploadedFiles['file'] = [(temp, prefix+'-processed'+suffix)]
    return showSimulationOptions()

def showSimulationOptions():
    return render_template('simulationOptions.html')

@app.route('/setSimulationOptions', methods=['POST'])
def setSimulationOptions():
    for key in request.form:
        session[key] = request.form[key]
    session['writeTrajectory'] = 'writeTrajectory' in request.form
    session['writeData'] = 'writeData' in request.form
    session['writeCheckpoint'] = 'writeCheckpoint' in request.form
    session['dataFields'] = request.form.getlist('dataFields')
    session['hmr'] = 'hmr' in request.form
    session['writeSimulationXml'] = 'writeSimulationXml' in request.form
    session['writeFinalState'] = 'writeFinalState' in request.form
    return createScript()

@app.route('/downloadScript')
def downloadScript():
    response = make_response(createScript())
    response.headers['Content-Disposition'] = 'attachment; filename="run_openmm_simulation.py"'
    return response

@app.route('/downloadPDB')
def downloadPDB():
    file, name = uploadedFiles['file'][0]
    file.seek(0, 0)
    response = make_response(file.read())
    response.headers['Content-Disposition'] = 'attachment; filename="%s"' % name
    return response

@app.route('/downloadPackage')
def downloadPackage():
    temp = tempfile.NamedTemporaryFile()
    with zipfile.ZipFile(temp, 'w', zipfile.ZIP_DEFLATED) as zip:
        zip.writestr('openmm_simulation/run_openmm_simulation.py', createScript())
        for key in uploadedFiles:
            for file, name in uploadedFiles[key]:
                file.seek(0, 0)
                zip.writestr('openmm_simulation/%s' % name, file.read())
    temp.seek(0, 0)
    return send_file(temp, 'application/zip', True, 'openmm_simulation.zip', max_age=0)

@app.route('/showRunSimulation')
def showRunSimulation():
    homeDir = os.path.expanduser('~')
    defaultDir = os.path.join(homeDir, 'openmm_simulation')
    return render_template('runSimulation.html', defaultDir=defaultDir)

@app.route('/startSimulation', methods=['POST'])
def startSimulation():
    global scriptOutput, simulationProcess
    conn1, conn2 = Pipe()
    scriptOutput = conn1
    # Create the simulation directory and copy files.
    try:
        outputDir = request.form['directory']
        if not os.path.isdir(outputDir):
            os.makedirs(outputDir)
    except:
        conn2.send('An error occurred while creating the simulation directory: %s' % sys.exc_info()[1])
        conn2.send(None)
        return ""
    try:
        for key in uploadedFiles:
            for file, name in uploadedFiles[key]:
                file.seek(0, 0)
                with open(os.path.join(outputDir, name), 'wb') as outputFile:
                    shutil.copyfileobj(file, outputFile)
        with open(os.path.join(outputDir, 'run_openmm_simulation.py'), 'w') as outputFile:
            outputFile.write(createScript())
    except:
        conn2.send('An error occurred while copying the input files: %s' % sys.exc_info()[1])
        conn2.send(None)
        return ""
    # Run the simulation in a subprocess.
    simulationProcess = Process(target=simulate, args=(conn2, outputDir, createScript(True)))
    simulationProcess.start()
    return ""

@app.route('/stopSimulation', methods=['POST'])
def stopSimulation():
    global scriptOutput, simulationProcess
    os.kill(simulationProcess.pid, signal.SIGKILL)
    scriptOutput = None
    return ""

@app.route('/getSimulationOutput')
def getSimulationOutput():
    global scriptOutput
    if scriptOutput is None:
        return "", 404
    output = []
    try:
        while scriptOutput.poll():
            data = scriptOutput.recv()
            if data is None:
                scriptOutput = None
                break
            else:
                output.append(data)
    except EOFError:
        scriptOutput = None
    return "".join(output)

def simulate(output, outputDir, script):
    try:
        exec(script, {"output":output, "outputDir":outputDir})
    except Exception as e:
        output.send('\nThe simulation failed with the following error:\n\n')
        output.send(str(e))
        output.send('\n\nDetails:\n\n')
        output.send(traceback.format_exc())
    output.send(None)

def configureDefaultOptions():
    """Select default options based on the file format and force field."""
    implicitWater = False
    if (session['fileType'] == 'pdb' or session['fileType'] == 'tinker') and session['waterModel'] == 'implicit':
        implicitWater = True
    isAmoeba = (session['fileType'] == 'pdb' and 'amoeba' in session['forcefield']) or (session['fileType'] == 'tinker')
    isDrude = session['fileType'] == 'pdb' and session['forcefield'].startswith('charmm_polar')
    session['isAmoeba'] = isAmoeba
    session['ensemble'] = 'nvt' if implicitWater else 'npt'
    session['platform'] = 'auto'
    session['precision'] = 'single'
    if isAmoeba:
        session['cutoff'] = '0.7'
        session['vdwCutoff'] = '0.9'
        session['inducedEpsilon'] = '0.00001'
    else:
        session['cutoff'] = '2.0' if implicitWater else '1.0'
    session['ewaldTol'] = '0.0005'
    session['constraintTol'] = '0.000001'
    session['hmr'] = True
    session['hmrMass'] = '1.5'
    if isAmoeba:
        session['dt'] = '0.002'
    elif isDrude:
        session['dt'] = '0.001'
    else:
        session['dt'] = '0.004'
    session['steps'] = '1000000'
    session['equilibrationSteps'] = '1000'
    session['temperature'] = '300'
    session['friction'] = '1.0'
    session['pressure'] = '1.0'
    session['barostatInterval'] = '25'
    if implicitWater:
        session['nonbondedMethod'] = 'NoCutoff' if isAmoeba else 'CutoffNonPeriodic'
    else:
        session['nonbondedMethod'] = 'PME'
    session['writeTrajectory'] = True
    session['trajFormat'] = 'dcd'
    session['trajFilename'] = 'trajectory.dcd'
    session['trajInterval'] = '10000'
    session['writeData'] = True
    session['dataFilename'] = 'log.txt'
    session['dataInterval'] = '1000'
    session['dataFields'] = ['step', 'speed' ,'progress', 'potentialEnergy', 'temperature']
    session['writeCheckpoint'] = True
    session['checkpointFilename'] = 'checkpoint.chk'
    session['checkpointInterval'] = '10000'
    session['writeSimulationXml'] = False
    session['systemXmlFilename'] = 'system.xml'
    session['integratorXmlFilename'] = 'integrator.xml'
    session['writeFinalState'] = False
    session['finalStateFileType'] = 'stateXML'
    session['finalStateFilename'] = "final_state.xml"
    if isAmoeba:
        session['constraints'] = 'none'
    else:
        session['constraints'] = 'hbonds'

def createScript(isInternal=False):
    script = []

    # If we are creating this script for internal use to run a simulation directly, add extra code at the top
    # to set the working directory and redirect stdout to the pipe.

    if isInternal:
        script.append("""
import os
import sys
import time

class PipeOutput(object):
    def write(self, string):
        output.send(string)

sys.stdout = PipeOutput()
sys.stderr = PipeOutput()
os.chdir(outputDir)""")

    # Header
    
    script.append('# This script was generated by OpenMM-Setup on %s.\n' % datetime.date.today())
    script.append('from openmm import *')
    script.append('from openmm.app import *')
    script.append('from openmm.unit import *')
    
    # Input files
    
    script.append('\n# Input Files\n')
    fileType = session['fileType']
    isAmoeba = session['isAmoeba']
    if fileType == 'pdb':
        pdbType = session['pdbType']
        if pdbType == 'pdb':
            script.append("pdb = PDBFile('%s')" % uploadedFiles['file'][0][1])
        else:
            script.append("pdbx = PDBxFile('%s')" % uploadedFiles['file'][0][1])
        forcefield = session['forcefield']
        water = session['waterModel']
        if forcefield == 'amoeba2018.xml':
            water = ('amoeba2018_gk.xml' if water == 'implicit' else None)
        elif forcefield == 'charmm_polar_2019.xml' or forcefield == 'charmm_polar_2023.xml':
            water = None
        elif water == 'implicit':
            models = {'amber99sb.xml': 'amber99_obc.xml',
                      'amber99sbildn.xml': 'amber99_obc.xml',
                      'amber03.xml': 'amber03_obc.xml',
                      'amber10.xml': 'amber10_obc.xml'}
            water = models[forcefield]
        if water is None:
            script.append("forcefield = ForceField('%s')" % forcefield)
        else:
            script.append("forcefield = ForceField('%s', '%s')" % (forcefield, water))
    elif fileType == 'amber':
        script.append("inpcrd = AmberInpcrdFile('%s')" % uploadedFiles['inpcrdFile'][0][1])
        script.append("prmtop = AmberPrmtopFile('%s', periodicBoxVectors=inpcrd.boxVectors)" % uploadedFiles['prmtopFile'][0][1])
    elif fileType == 'charmm':
        script.append("psf = CharmmPsfFile('%s')" % uploadedFiles['psfFile'][0][1])
        script.append("crd = CharmmCrdFile('%s')" % uploadedFiles['crdFile'][0][1])
        ffFiles = ', '.join(["'%s'" % f[1] for f in uploadedFiles['ffFiles']])
        script.append("params = CharmmParameterSet(%s)" % ffFiles)
    elif fileType == 'gromacs':
        script.append("gro = GromacsGroFile('%s')" % uploadedFiles['groFile'][0][1])
        script.append("top = GromacsTopFile('%s', includeDir='%s'," % (uploadedFiles['topFile'][0][1], session['gromacsIncludeDir']))
        script.append("    periodicBoxVectors=gro.getPeriodicBoxVectors())")
    elif fileType == 'tinker':
        ffFiles = ', '.join(["'%s'" % f[1] for f in uploadedFiles['ffFiles']])
        script.append("tinker = TinkerFiles('%s', [%s])" % (uploadedFiles['xyzFile'][0][1], ffFiles))

    # System configuration

    script.append('\n# System Configuration\n')
    nonbondedMethod = session['nonbondedMethod']
    script.append('nonbondedMethod = %s' % nonbondedMethod)
    if isAmoeba:
        if nonbondedMethod != 'NoCutoff':
            script.append('multipoleCutoff = %s*nanometers' % session['cutoff'])
            script.append('vdwCutoff = %s*nanometers' % session['vdwCutoff'])
        script.append('mutualInducedTargetEpsilon = %s' % session['inducedEpsilon'])
        constraints = 'none'
    else:
        if nonbondedMethod != 'NoCutoff':
            script.append('nonbondedCutoff = %s*nanometers' % session['cutoff'])
        constraints = session['constraints']
        constraintMethods = {'none': 'None',
                             'water': 'None',
                             'hbonds': 'HBonds',
                             'allbonds': 'AllBonds'}
        script.append('constraints = %s' % constraintMethods[constraints])
        script.append('rigidWater = %s' % ('False' if constraints == 'none' else 'True'))
        if constraints != 'none':
            script.append('constraintTolerance = %s' % session['constraintTol'])
        if session['hmr']:
            script.append('hydrogenMass = %s*amu' % session['hmrMass'])
    if nonbondedMethod == 'PME':
        script.append('ewaldErrorTolerance = %s' % session['ewaldTol'])

    # Integration options

    script.append('\n# Integration Options\n')
    script.append('dt = %s*picoseconds' % session['dt'])
    script.append('temperature = %s*kelvin' % session['temperature'])
    script.append('friction = %s/picosecond' % session['friction'])
    ensemble = session['ensemble']
    if ensemble == 'npt':
        script.append('pressure = %s*atmospheres' % session['pressure'])
        script.append('barostatInterval = %s' % session['barostatInterval'])

    # Simulation options

    script.append('\n# Simulation Options\n')
    script.append('steps = %s' % session['steps'])
    script.append('equilibrationSteps = %s' % session['equilibrationSteps'])
    if session['platform'] != 'auto':
        script.append("platform = Platform.getPlatformByName('%s')" % session['platform'])
        if session['platform'] in ('CUDA', 'HIP', 'OpenCL'):
            script.append("platformProperties = {'Precision': '%s'}" % session['precision'])
    if session['writeTrajectory']:
        if session['trajFormat'] == 'dcd':
            script.append("dcdReporter = DCDReporter('%s', %s)" % (session['trajFilename'], session['trajInterval']))
        else:
            script.append("xtcReporter = XTCReporter('%s', %s)" % (session['trajFilename'], session['trajInterval']))
    if session['writeData']:
        args = ', '.join('%s=True' % field for field in session['dataFields'])
        script.append("dataReporter = StateDataReporter('%s', %s, totalSteps=steps," % (session['dataFilename'], session['dataInterval']))
        script.append("    %s, separator='\\t')" % args)
        if isInternal:
            # Create a second reporting sending to stdout so we can display it in the browser.
            script.append("consoleReporter = StateDataReporter(sys.stdout, %s, totalSteps=steps, %s, separator='\\t')" % (session['dataInterval'], args))
    if session['writeCheckpoint']:
        script.append("checkpointReporter = CheckpointReporter('%s', %s)" % (session['checkpointFilename'], session['checkpointInterval']))
    
    # Prepare the simulation
    
    script.append('\n# Prepare the Simulation\n')
    script.append("print('Building system...')")
    if fileType == 'pdb' and pdbType == 'pdb':
        script.append('topology = pdb.topology')
        script.append('positions = pdb.positions')
    elif fileType == 'pdb' and pdbType == 'pdbx':
        script.append('topology = pdbx.topology')
        script.append('positions = pdbx.positions')
    elif fileType == 'amber':
        script.append('topology = prmtop.topology')
        script.append('positions = inpcrd.positions')
    elif fileType == 'charmm':
        script.append('topology = psf.topology')
        script.append('positions = crd.positions')
    elif fileType == 'gromacs':
        script.append('topology = top.topology')
        script.append('positions = gro.positions')
    elif fileType == 'tinker':
        script.append('topology = tinker.topology')
        script.append('positions = tinker.positions')
    if fileType == 'pdb' and not isAmoeba:
        if forcefield == 'charmm_polar_2019.xml' or forcefield == 'charmm_polar_2023.xml' or 'tip4p' in water or 'tip5p' in water:
            script.append('modeller = Modeller(topology, positions)')
            script.append('modeller.addExtraParticles(forcefield)')
            script.append('topology = modeller.topology')
            script.append('positions = modeller.positions')
    hmrOptions = ', hydrogenMass=hydrogenMass' if session['hmr'] else ''
    if isAmoeba:
        cutoffOptions = ' nonbondedCutoff=multipoleCutoff, vdwCutoff=vdwCutoff,' if nonbondedMethod == 'PME' else ''
    else:
        cutoffOptions = ' nonbondedCutoff=nonbondedCutoff,' if nonbondedMethod != 'NoCutoff' else ''
    if fileType  == 'pdb':
        script.append('system = forcefield.createSystem(topology, nonbondedMethod=nonbondedMethod,%s' % cutoffOptions)
        if isAmoeba:
            script.append('    mutualInducedTargetEpsilon=mutualInducedTargetEpsilon)')
        else:
            script.append('    constraints=constraints, rigidWater=rigidWater%s%s)' % (', ewaldErrorTolerance=ewaldErrorTolerance' if nonbondedMethod == 'PME' else '', hmrOptions))
    elif fileType == 'amber':
        script.append('system = prmtop.createSystem(nonbondedMethod=nonbondedMethod,%s' % cutoffOptions)
        script.append('    constraints=constraints, rigidWater=rigidWater%s%s)' % (', ewaldErrorTolerance=ewaldErrorTolerance' if nonbondedMethod == 'PME' else '', hmrOptions))
    elif fileType == 'charmm':
        script.append('system = psf.createSystem(params, nonbondedMethod=nonbondedMethod,%s' % cutoffOptions)
        script.append('    constraints=constraints, rigidWater=rigidWater%s%s)' % (', ewaldErrorTolerance=ewaldErrorTolerance' if nonbondedMethod == 'PME' else '', hmrOptions))
    elif fileType == 'gromacs':
        script.append('system = top.createSystem(nonbondedMethod=nonbondedMethod,%s' % cutoffOptions)
        script.append('    constraints=constraints, rigidWater=rigidWater%s%s)' % (', ewaldErrorTolerance=ewaldErrorTolerance' if nonbondedMethod == 'PME' else '', hmrOptions))
    elif fileType == 'tinker':
        cutoffOptions = ' nonbondedCutoff=multipoleCutoff, vdwCutoff=vdwCutoff,' if nonbondedMethod == 'PME' else ''
        ewaldOptions = ', ewaldErrorTolerance=ewaldErrorTolerance' if nonbondedMethod == 'PME' else ''
        waterOptions = ', implicitSolvent=True' if session['waterModel'] == 'implicit' else ''
        script.append('system = tinker.createSystem(nonbondedMethod=nonbondedMethod,%s' % cutoffOptions)
        script.append('    mutualInducedTargetEpsilon=mutualInducedTargetEpsilon%s%s)' % (ewaldOptions, waterOptions))
    if ensemble == 'npt':
        script.append('system.addForce(MonteCarloBarostat(pressure, temperature, barostatInterval))')
    if isAmoeba:
        # Use a MTSLangevinIntegrator.
        script.append('for force in system.getForces():')
        script.append('    if isinstance(force, AmoebaMultipoleForce) or isinstance(force, AmoebaVdwForce) or isinstance(force, AmoebaGeneralizedKirkwoodForce):')
        script.append('        force.setForceGroup(1)')
        script.append('integrator = MTSLangevinIntegrator(temperature, friction, dt, [(0,2), (1,1)])')
    elif fileType == 'pdb' and forcefield.startswith('charmm_polar'):
        script.append('integrator = DrudeLangevinIntegrator(temperature, friction, 1*kelvin, friction, dt)')
    else:
        script.append('integrator = LangevinMiddleIntegrator(temperature, friction, dt)')
    if constraints != 'none':
        script.append('integrator.setConstraintTolerance(constraintTolerance)')
    script.append('simulation = Simulation(topology, system, integrator%s%s)' % (', platform' if session['platform'] != 'auto' else '', ', platformProperties' if session['platform'] in ('CUDA', 'HIP', 'OpenCL') else ''))
    script.append('simulation.context.setPositions(positions)')

    # Output XML files for system and integrator

    if session['writeSimulationXml']:
        def _xml_script_segment(to_serialize, target_file):
            if target_file == "":
                # if filename is blank, we cannot create the file
                return []
            return [
                f'with open("{target_file}", mode="w") as file:',
                f'    file.write(XmlSerializer.serialize({to_serialize}))'
            ]

        script.append("\n# Write XML serialized objects\n")
        script.extend(_xml_script_segment('system', session['systemXmlFilename']))
        script.extend(_xml_script_segment('integrator', session['integratorXmlFilename']))

    # Minimize and equilibrate
    
    script.append('\n# Minimize and Equilibrate\n')
    script.append("print('Performing energy minimization...')")
    script.append('simulation.minimizeEnergy()')
    script.append("print('Equilibrating...')")
    script.append('simulation.context.setVelocitiesToTemperature(temperature)')
    script.append('simulation.step(equilibrationSteps)')
    
    # Simulate
    
    script.append('\n# Simulate\n')
    script.append("print('Simulating...')")
    if session['writeTrajectory']:
        if session['trajFormat'] == 'dcd':
            script.append('simulation.reporters.append(dcdReporter)')
        else:
            script.append('simulation.reporters.append(xtcReporter)')
    if session['writeData']:
        script.append('simulation.reporters.append(dataReporter)')
        if isInternal:
            script.append('simulation.reporters.append(consoleReporter)')
    if session['writeCheckpoint']:
        script.append('simulation.reporters.append(checkpointReporter)')
    script.append('simulation.currentStep = 0')
    script.append('simulation.step(steps)')

    # Output final simulation state
    if session['writeFinalState']:
        script.append("\n# Write file with final simulation state\n")
        state_script = {
            'checkpoint': ['simulation.saveCheckpoint("{filename}")'],
            'stateXML': ['simulation.saveState("{filename}")'],
            'pdbx': ['state = simulation.context.getState(getPositions=True, enforcePeriodicBox=system.usesPeriodicBoundaryConditions())',
                     'with open("{filename}", mode="w") as file:',
                     '    PDBxFile.writeFile(simulation.topology, state.getPositions(), file)'],
        }[session['finalStateFileType']]
        lines = [line.format(filename=session['finalStateFilename']) for line in state_script]
        script.extend(lines)

    return "\n".join(script)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-p', '--port', type=int, default=5000)
    port = parser.parse_args().port
    if port < 1 or port > 65535:
        raise ValueError(f"The given {port=} is not between 1 and 65535")

    def open_browser():
        # Give the server a moment to start before opening the browser.
        time.sleep(1)
        url = f'http://127.0.0.1:{port}'
        webbrowser.open(url)

    global server, shutdownEvent
    server = make_server('localhost', port, app)
    shutdownEvent = threading.Event()
    threading.Thread(target=open_browser).start()
    threading.Thread(target=server.serve_forever).start()
    shutdownEvent.wait()
    server.shutdown()

if __name__ == '__main__':
    main()
