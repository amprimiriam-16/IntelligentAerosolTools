########################################################################### 
# 
#  'SMILES_to_AIOMFAC_input.py'
#
#  :: Purpose ::
#  Python script to process a list of organic molecular structures in the 
#  form of SMILES read from an input file to determine the subgroup 
#  description of the molecules for use as input with the AIOMFAC model.
#  The determined subgroup fingerprints of the molecules are stored in array
#  cpsubs and an AIOMFAC-web style input file is written based on this.
#  By default, the AIOMFAC input file will contain also water (see param.).
#
#  :: Requirements ::
#  Aside from a python installation, the epam indigo toolkit needs to be 
#  installed, e.g. via command line:  pip install epam.Indigo
#
#  :: Specific Settings ::
#  Some specific parameters, e.g. for the treatment of radical species, can 
#  be set near the top of this file in the Switches & Parameters section.
#
#  :: Authors & Copyright ::                                         
#  Dalrin Ampritta Amaladhasan, Andreas Zuend, 
#  Dept. Atmospheric and Oceanic Sciences, McGill University                                    
#
###########################################################################

from indigo import *
from IndigoToolkit.indigo_renderer import *
import ModSmilesTools
from SMARTS_query_list import SMARTS_AIOMFAC
import sys
import time
#import argparse

#--- Switches & Parameters ----------------------------------------------------
debugging_verbose = False   #set to True to write detailed information to the terminal window and creation of images; set False to avoid most output;
replaceRadicals = True     #if set True, component SMILES indicating radical atoms (e.g. -C[O] or -C[O.]) will be 
                           #modified into non-radical species (by adding implicit H to reach non-radical valence);
replaceRdbOO = True        #if set True, replace SMILES containing C=[O+][O-] with a similar one with C-O-OH hydroperoxy acid group;
waterAsComp01 = True       #add water to the generated AIOMFAC input file (while H2O is not part of the SMILES list);
#------------------------------------------------------------------------------

#initialize some global variables:
timeStart = time.perf_counter()
indigo = Indigo()
renderer = IndigoRenderer(indigo)
SMILESlist = []
gpfile2 = []

print(' ')
cmdline = sys.argv[0:]
if len(cmdline) > 1:
    smiles_file = sys.argv[1]
    debugging_verbose = False
else: 
    # for debugging, use the following input file path instead of the command line input
    #smiles_file = "InputFiles/smiles_1250.txt"  #"InputFiles/smiles_1441.txt"
    smiles_file = "InputFiles/smiles_1415.txt"

if debugging_verbose:
   debugRenderedStructures = True   #set to True for debugging and showing of highlighted subgroup in a molecule (file will be overwritten);
   showRenderedStructures = True    #set to True for debugging and/or graphical output of png images of the structures using the indigo renderer;
else:
   debugRenderedStructures = False
   showRenderedStructures = False 

if debugging_verbose:
   print('The relative path and name of the SMILES input file is : ', str(smiles_file))
   print('')

#determine the relative path to output folder using the location of this  *.py file as the starting point.
#use of the os.path functions is necessary to ensure proper path strings when calling this .py file 
#from a non-local directory, e.g. the 2D-lumping framework program;
pyfilepath = os.path.abspath(cmdline[0])
locpath = os.path.dirname(pyfilepath)
outpath = os.path.abspath(locpath +'./OutputFiles')

#clean up output figure folder to avoid presence of misleading older files:
dirfilelist = os.listdir(outpath)
for structfile in dirfilelist:
   if structfile.endswith('.png') and structfile.index('Structure_') != None:
      os.remove(os.path.join(outpath, structfile))

AbsPathFName = os.path.abspath(smiles_file)
with open(AbsPathFName, "r") as inpfile:
   for line in inpfile:
      SMILESlist.append(line)
 
#initialize lists:
ii = len(SMILESlist)
is_invalidSMILES = [True]*ii
is_PureAlcohol = [False]*ii
istatus = [-1]*ii
highLAtomsList = []        #list of highlighted atoms
igatoms_index = []         #list of ignored atoms in a structure (used in the indigo substructure matching below)
cpsubs = [[0 for col in range(0,173)] for row in range(len(SMILESlist))]   #2D array with as many rows as rows in SMILESlist and 172 columns (for the 172 different AIOMFAC subgroups) 

# alkylPlusOH = list with the CHn groups attached to hydroxyl:
alkylOH = ['[CH3][OX2H1]', '[CH2][OX2H1]', '[CH1][OX2H1]', '[CH0][OX2H1]']
subgCHnOH = [149, 150, 151, 152]    #(there is also AIOMFAC subgroup 153 = -OH present)

# alkylTailEnd = list of SMARTS for detecting the two end carbons of a hydrophobic tail in alcohols:
alkylTailEnd = ['[CH3][CH2]', '[CH3][CH1]', '[CH3][CH0]']
subgTailEnd  = [146, 147, 148]      #each of these matches will also include subgroup no. 145 for -CH3; 

# alkylGroups = list of SMARTS for general alkyl subgroups for use with alcohols:
alkylG = ['[CH3]', '[CH2]', '[CH1]', '[CH0]']
subgTail = [145, 146, 147, 148]     #list of special alkyl in a hydrophobic tail AIOMFAC subgroup numbers
subgAlc  = [141, 142, 143, 144]     #list of other "in alcohol" alkyl AIOMFAC subgroup numbers;

#====== Function ===========================================================
#Define a function that works on the main task of matching SMARTS patterns and returning the subgroup array;
def detAIOMFACsubgs(ind, molec, is_PureAlc):
   subgs = [0 for j in range(0,173)] 
   igatoms_index.clear()  
   structure = indigo.loadMolecule(molec) 
   if structure.canonicalSmiles() != '':
      matcher = indigo.substructureMatcher(structure)
   else:
      info = -99
      return subgs, info   #exceptional return point;

   if debugging_verbose:
      print('')
      print('Input SMILES of molecule: ' + structure.canonicalSmiles() + ' ind ' , ind)
   molecatomcount = structure.countAtoms()

   if (not is_PureAlc):
      #======
      # (1) iterate over all general AIOMFAC SMARTS patterns and count number of matches in current molecule:
      alkyldeduct = 0
      CH3deduct = 0
      for row in range(len(SMARTS_AIOMFAC[:])):
         sub,item = SMARTS_AIOMFAC[row][0:]
         tsmarts = indigo.loadSmarts(item)    
         match1 = matcher.match(tsmarts)
         has_SMARTS_match = False
         while match1 != None: #then match found... 
            #...check for exception cases and count matched subgroup:
            if sub == 1550:
               subgs[155] = subgs[155] +1
               alkyldeduct += 1     #counter for how many alkyl groups (detected later) should be deducted due to an exception case match here;
            elif sub == 1580:
               subgs[158] = subgs[158] +1
               alkyldeduct += 1
            elif sub == 1700:
               subgs[170] = subgs[170] +1
               alkyldeduct += 1
            elif sub == 1710:
               subgs[171] = subgs[171] +1
               alkyldeduct += 2 
            elif sub == 2002:
               subgs[20] = subgs[20] +1
               alkyldeduct += 1 
            elif sub == 2022:
               subgs[20] = subgs[20] +1
               subgs[22] = subgs[22] +1
            elif sub == 2220:
               subgs[20] = subgs[20] +1
               subgs[22] = subgs[22] +1
               alkyldeduct += 1
            elif sub == 2224:
               subgs[20] = subgs[20] +1
               subgs[24] = subgs[24] +1
               CH3deduct += 1
            elif sub == 2502:
               subgs[25] = subgs[25] +1
               alkyldeduct += 1
            else: #regular case
               subgs[sub] = subgs[sub] +1
            #...
            for atom in tsmarts.iterateAtoms():
               if match1.mapAtom(atom) != None:
                  atom_index = match1.mapAtom(atom).index()
                  igatoms_index.append(atom_index)
                  atom_to_ignore = structure.getAtom(atom_index)
                  if debugRenderedStructures:
                     has_SMARTS_match = True
                     atom_to_ignore.highlight()          #just for debugging plots; remove otherwise
                  matcher.ignoreAtom(atom_to_ignore) 
            match1 = matcher.match(tsmarts)              #load the matcher again since some atoms from this molecule may have been ignored;

         #---- save an png-image of this component with highlighted subgroup (for debugging, otherwise not needed):
         if debugRenderedStructures and has_SMARTS_match:
            #produce output file with highlighted subgroup matches, all else in black:
            indigo.setOption("render-output-format", "png")
            indigo.setOption("render-margins", "5, 5")
            indigo.setOption("render-coloring", "False")
            indigo.setOption("render-background-color", "1.0, 1.0, 1.0")
            indigo.setOption("render-atom-ids-visible", "True")
            indigo.setOption("render-label-mode", "Hetero")
            indigo.setOption("render-highlight-color-enabled", "True")
            indigo.setOption("render-highlight-color", "1.0, 0.65, 0.18")    #lighter orange for highlighted atoms
            st = str(ind).zfill(3)
            csub = str(sub).zfill(3)
            crow = str(row).zfill(2)  #to mark the SMARTS priority position
            outpathfile = "OutputFiles/Structure_" +st +"_SMARTS_" +crow +"_subgr_" +csub +".png"
            renderer.renderToFile(structure, outpathfile)
            #now unhighlight the atoms from the current match:
            for atom in structure.iterateAtoms():
               atom.unhighlight()

      #check whether CH3 alkyl groups should be deducted (e.g. due to imperfect matching of perester group):
      if CH3deduct > 0:
         deduct_is_possible = True
         while CH3deduct > 0 and deduct_is_possible:
            if subgs[1] > 0:
               subgs[1] -= 1
               CH3deduct -= 1
            else:
               deduct_is_possible = False
               alkyldeduct += CH3deduct
      #check whether any alkyl groups should be deducted because of imperfect subgroup matching by other groups;
      #(doing this will ensure that the number of C atoms is at least correctly accounted for):
      if alkyldeduct > 0:
         deduct_is_possible = True
         while alkyldeduct > 0 and deduct_is_possible:
            if subgs[3] > 0:
               subgs[3] -= 1
               alkyldeduct -= 1
            elif subgs[2] > 0:
               subgs[2] -= 1
               alkyldeduct -= 1
            elif subgs[4] > 0:
               subgs[4] -= 1
               alkyldeduct -= 1
            elif subgs[1] > 0:
               subgs[1] -= 1
               alkyldeduct -= 1
            else:
               deduct_is_possible = False
      #=======
   else:  #--> it is a pure aliphatic alcohol/polyol, so perform special pattern searches (steps 1 - 3);
      #--------
      # step (1) determine the CHn with OH groups and the OH groups themselves;
      for i,item in enumerate(alkylOH): 
         tsmarts = indigo.loadSmarts(item)    
         matchCHnOH = matcher.match(tsmarts)
         while matchCHnOH != None:
            sub = subgCHnOH[i]
            subgs[sub] = subgs[sub] +1   #count the matched CHn (with OH) subgroup; 
            subgs[153] = subgs[153] +1   #also count the matched OH subgroup;
            #...
            for atom in tsmarts.iterateAtoms():
               atom_to_ignore = matchCHnOH.mapAtom(atom)
               igatoms_index.append(atom_to_ignore.index())
               matcher.ignoreAtom(atom_to_ignore)
            #...
            matchCHnOH = matcher.match(tsmarts)          #load the matcher again since some atoms from this molecule may have been ignored;
      #--------
      # step (2) determine hydrophobic chains that terminate in -CHn-CH3  (n = 0,1,2); 
      #    [maybe use Indigo highlight and isHighlighted functions to mark CHn atoms that neighbor tail groups for if else testing in code to see whether chain of hydrophobic tail continues.]

      highLAtomsList.clear()
      # (2a) determine the hydrophobic tail end groups (if present) and highlight neighboring atoms (used later as tail marker);
      for i,item in enumerate(alkylTailEnd): 
         tsmarts = indigo.loadSmarts(item)    
         matchCHntail = matcher.match(tsmarts)
         while matchCHntail != None:
            sub = subgTailEnd[i]
            subgs[sub] = subgs[sub] +1      #register the matched CHn (n = 0,1,2) subgroup; 
            subgs[145] = subgs[145] +1      #also register the matched CH3 in hydrophobic tail subgroup;
            #...
            for atom in tsmarts.iterateAtoms():
               targetAtom = matchCHntail.mapAtom(atom)
               # first highlight the neighboring atom of the matched group, then ignore the matched atoms:
               igatoms_index.append(targetAtom.index())
               for neiatom in targetAtom.iterateNeighbors():
                  neiatom.highlight()
                  neiatom_index = neiatom.index()
               matcher.ignoreAtom(targetAtom)
            #...
            matchCHntail = matcher.match(tsmarts)           #load the matcher again since some atoms from this molecule may have been ignored;

      # check on the indices of highlighted atoms in the molecule:
      for atom in structure.iterateAtoms():
         if atom.index() in igatoms_index: 
            atom.unhighlight()   #exclude ignored atoms from being highlighted (here only show highlighted neighbour atoms that have not been mapped to subgroups);
         elif atom.isHighlighted():
            highLAtomsList.append(atom.index())

      # (2b) highlight all other hydrophobic tail atoms and associate them with specific alkyl tail subgroups:
      temporaryIgnoreL = []
      tsmarts = indigo.loadSmarts('[CX4]')                  #this SMARTS will hit any alkyl carbon that has not been ignored;
      matchCHntail = matcher.match(tsmarts)
      while matchCHntail != None:
         #...
         for atom in tsmarts.iterateAtoms():
            targetAtom = matchCHntail.mapAtom(atom)
            if targetAtom.isHighlighted():                  #select highlighted atoms only!
               #-.-
               for neiatom in targetAtom.iterateNeighbors():
                  neiatom_index = neiatom.index()
                  if neiatom_index not in igatoms_index and (not neiatom.isHighlighted()):
                     neiatom.highlight()
                     highLAtomsList.append(neiatom_index)
                     if neiatom_index in temporaryIgnoreL:
                        temporaryIgnoreL.remove(neiatom_index)
                        matcher.unignoreAtom(neiatom)       #unignore this atom as it is now highlighted and a tail atom;
               #-.-
               #determine the specific tail alkyl subgroups and then ignore the atom:
               kH = targetAtom.countHydrogens()
               sub = 148 - kH                               #here to hit the correct subgroup ID for hydrophobic tail CHn groups
               subgs[sub] = subgs[sub] +1
               igatoms_index.append(targetAtom.index())
               matcher.ignoreAtom(targetAtom)               #ignore current atom to allow matcher to progress to next hit;
            else:
               temporaryIgnoreL.append(targetAtom.index())  #add atom index of matched atom to a temporary ignore list; later unignored for final matching;
               matcher.ignoreAtom(targetAtom)               #ignore current atom to allow matcher to progress to next hit;
         #...
         matchCHntail = matcher.match(tsmarts)              #load the matcher again since some atoms from this molecule may have been ignored;
      
      #unignore temporarily ignored atoms:
      for atom in structure.iterateAtoms():
         if atom.index() in temporaryIgnoreL: 
            matcher.unignoreAtom(atom) 
            temporaryIgnoreL.remove(atom.index())

      #--------
      # step (3) assign alkyl within alcohols type to all remaining CHn groups.
      for i,item in enumerate(alkylG): 
         tsmarts = indigo.loadSmarts(item)    
         matchCHnAlc = matcher.match(tsmarts)
         while matchCHnAlc != None:
            sub = subgAlc[i]
            subgs[sub] = subgs[sub] +1   #register the matched CHn (n = 0,1,2) subgroup; 
            #...
            for atom in tsmarts.iterateAtoms():
               targetAtom = matchCHnAlc.mapAtom(atom)
               igatoms_index.append(targetAtom.index())
               matcher.ignoreAtom(targetAtom)
            #...
            matchCHnAlc = matcher.match(tsmarts)         #load the matcher again since some atoms from this molecule may have been ignored;

   #---- save an png-image of this component (for debugging)
   if showRenderedStructures:
      #produce new output file:
      indigo.setOption("render-output-format", "png")
      indigo.setOption("render-margins", "5, 5")
      indigo.setOption("render-coloring", "True")
      indigo.setOption("render-background-color", "1.0, 1.0, 1.0")
      indigo.setOption("render-atom-ids-visible", "False")
      indigo.setOption("render-label-mode", "Hetero")
      indigo.setOption("render-highlight-color-enabled", "True")
      indigo.setOption("render-highlight-color", "1.0, 0.65, 0.18")    #lighter orange for highlighted atoms
      st = str(ind).zfill(3)
      outpathfile = "OutputFiles/Structure_" +st +".png"
      renderer.renderToFile(structure, outpathfile)      #in this case, the highlighted orange atoms are those from the hydrophobic tail search  
                                                         #in step (2b), but not including the CH3-CHn- end groups.

   if debugging_verbose:
      print('Total number of atoms in input molecule are : ', molecatomcount)
      print('Total number of atoms that have been ignored are : ', len(igatoms_index))
      if subgs[0] > 0:
         print('WARNING: not all atoms were assigned to proper AIOMFAC subgroups!')
         print('The count of subgroup-0 atoms is non-zero (though it should always be 0).')
         input()
   #check whether issues occurred that may need attention:
   ii = molecatomcount - len(igatoms_index)
   if ii != 0:
      print('WARNING: not all atoms were correctly matched to AIOMFAC subgroups for the input SMILES at ind ' + str(ind).zfill(2))
      print('Number of non-H atoms not matched to any AIOMFAC subgroup: ' + str(ii).zfill(2))
      i = ii
   return subgs, ii
#=== end of function ===================================================================


#+++ main task: processing list of SMILES +++++++++++++++++++++++++++++++++++++++++++
for ind,molecule in enumerate(SMILESlist):
   st = str(molecule)
   #(1) run string processing on input SMILES to remove unneccessary attributes like / or \ or revise radical atoms; 
   # they are unresolved for AIOMFAC description of structures.
   st = st.replace('/','')
   st = st.replace('\\','')
   if replaceRadicals:           #replacing only radical atoms, but not charged atoms;
      st = st.replace('[O]','O')
      st = st.replace('[O.]','O')
      st = st.replace('[C]','C')
      st = st.replace('[C.]','C')
      st = st.replace('[N]','N')
      st = st.replace('[N.]','N')
      st = st.replace('[S]','S')
      st = st.replace('[S.]','S')
   if replaceRdbOO:
      st = st.replace('=[O+][O-]','OO')
      st = st.replace('[O-][O+]=','OO') 
   SMILESlist[ind] = st.strip()    #strip() removes leading and trailing whitespace, tabs, etc.
   molec = SMILESlist[ind]

   #(2) use special smarts to check validity of SMILES and whether a molecule contains 
   #    a hydroxyl group and no other/different non-carbon functionalities; 
   validS, errortext = ModSmilesTools.verifySMILES(molec)
   is_invalidSMILES[ind] = not validS
   is_PureAlcohol[ind] = ModSmilesTools.PureAlcoholCheck(molec)
   info = -77
   if (is_invalidSMILES[ind]):
      print('')
      print('input SMILES is: ' + indigo.loadMolecule(molec).canonicalSmiles() )
      print('WARNING: issues found with SMILES: ' + errortext)    #for debugging
      print('continue anyways... or manually stop execution and correct input file') 
      input()
   elif len(SMILESlist[ind]) < 1:
      is_invalidSMILES[ind] = True
   else:
      #(3) call function for main task of matching AIOMFAC subgroups via SMARTS:
      cpsubs[ind][:], istatus[ind] = detAIOMFACsubgs(ind, molec, is_PureAlcohol[ind])

   if ind > 1 and ind % 1000 == 0:  #modulo test passed;
      print(f'status update: just processed SMILES at index: {ind}')

#clean up output array in case of gaps in the SMILES list:
if any(is_invalidSMILES):
   ind = is_invalidSMILES.index(True)
else:
   ind = -1
while ind > -1:
   del is_invalidSMILES[ind]
   del SMILESlist[ind]           #remove item/row in list
   del cpsubs[ind]
   del is_PureAlcohol[ind]
   del istatus[ind]
   if any(is_invalidSMILES):
      ind = is_invalidSMILES.index(True)
   else:
      ind = -1
#+++++ end of main task +++++++++++++++++++++++++++++++++++++++++++

#save original SMILES input file with extension '_orig.txt' if it does not exist yet 
#and the modified SMILES list as a smiles_XXXX.txt file:
from shutil import copy
AbsPathFNameOrig = AbsPathFName.replace('.txt','_orig.txt')
if os.path.isfile(AbsPathFNameOrig) != True:
   copy(AbsPathFName, AbsPathFNameOrig)
fileA = open(AbsPathFName, 'w')
for item in SMILESlist:
   fileA.write("{}\n".format(item))
fileA.close()

#generate an 'AIOMFAC input file' name in folder ./OutputFiles:
head, smilesfilename1 = os.path.split(smiles_file)
AIOMFACfname = smilesfilename1.replace('smiles','input')
# use function writeAIOMFACfile with above data:
_ = ModSmilesTools.writeAIOMFACfile(AIOMFACfname, outpath, SMILESlist, cpsubs, waterAsComp01)

timeStop = time.perf_counter()

print('')
print('SMILES processing for AIOMFAC done. The generated AIOMFAC input file,', AIOMFACfname, ', is located')
print(f'in directory: {outpath}')
print(f'The SMILES processing and file creation took a time of ~ {timeStop - timeStart:0.1f} sec.')
print('')