#!/usr/bin/env python3.7
"""
Author: Mikhail K.
Pacific Lutheran University NSSURP 2020
Department of Computer Science

Usage:
python main.py -i ./input_folder_name -o output_folder_name

The input folder must contain at least one protein model (PDB format) for this tool to work.
"""

import os, sys, math, time, pwd
import argparse, subprocess
import numpy as np
from datetime import timedelta
from collections import namedtuple as NamedTuple
from os.path import join
from joblib import load
start_time = time.monotonic()

################################################################################
# GLOBALS ----------------------------------------------------------------------
###############################G#################################################
QA_DIR = "/path/to/SynthQA" # the main directory of this tool
ANDIS_PATH = "/path/to/ANDIS"
NUM_FEATURES = 14 # the number of non-TopQA features we trained on, do not modify unless you're training your own models.

################################################################################
# ARGUMENT PARSING -------------------------------------------------------------
################################################################################
parser = argparse.ArgumentParser()

parser._optionals.title = "keyword arguments"
parser.add_argument("-i", dest="input_dir", default = '', help = "input directory containing PDBs to score")
parser.add_argument("-o", dest="output_dir", default = '', help = "name of output directory folder")

if len(sys.argv) < 4:
    print("\n#-------- [NOT ENOUGH ARGUMENTS] --------#")
    parser.print_help(sys.stderr)
    print("#---------------------------------------#")
    sys.exit(1)

options = parser.parse_args()

# Check input PDB directory
if not (os.path.isdir(options.input_dir) and len(os.listdir(options.input_dir)) > 0):
    print("ERROR: Please check your input directory! It is either empty or doesn't exist.")
    sys.exit()

# Set argument variables
PDB_DIR = os.path.abspath(options.input_dir)
OUTPUT_DIR = os.path.abspath(options.output_dir)
CWD = os.getcwd()





################################################################################
# HELPER FUNCTIONS -------------------------------------------------------------
################################################################################

def assert_directory(dir):
    """ Makes a new directory if the one given doesn't exist. """
    if not os.path.isdir(dir):
        os.makedirs(dir)
    return dir

def sigmoid(x):
    """ Normalizes X using the sigmoid function. """
    if x < 0:
        return 1 - 1/(1 + math.exp(x))
    else:
        return 1/(1 + math.exp(-x))

def rescale(x, min_bound, max_bound):
  """
  Normalizes X by scaling it to a value between 0 and 1 that is
  proportionally relevant to min_bound and max_bound.

  min_bound and max_bound correspond to the minimum and maximum values of the unnormalized
  feature in our training data.
  """
  if x < min_bound:
    x = min_bound
  elif x > max_bound:
    x = max_bound
  else:
    x = x
  new_value = (x - min_bound) * (1/(max_bound-min_bound))
  assert (new_value <= 1) and (new_value >= 0), "Something went wrong. You shouldn't be seeing this message!" # make sure that it's between 1 and 0
  return new_value



################################################################################
# FEATURES ---------------------------------------------------------------------
################################################################################
# Each function in this section is responsible for executing scripts           #
# that generate features, which are saved in their own respective files in the #
# output directory. The format is either a new file for each pdb               #
# (e.g server1.euc, server2.euc, server3.euc...), or a single output file where#
# each row represents a single pdb (e.g. surface_score.txt)                    #
################################################################################

# TODO: replace every filepath-related string with the os.path function equivelant (e.g. join() chdir(), etc.)
def gen_euclidean():
    """
    Generates a Euclidean Compact Score for each input file.
    Results for each model are saved as 'OUTPUT_DIR/euclidean/MODEL_NAME.euc'
    """
    output_path = assert_directory(OUTPUT_DIR + "/euclidean")
    log_folder = assert_directory(OUTPUT_DIR + "/log")
    for pdbfile in os.listdir(PDB_DIR):
        outfile = join(output_path, pdbfile.rstrip(".pdb")+'.euc')
        pdbfile = join(PDB_DIR, pdbfile)
        script = join(QA_DIR, "scripts/Euclidean_compact.pl")
        f = open(log_folder+"/euclidean.log", 'a+')
        euc_proc = subprocess.Popen(["perl", script, pdbfile, outfile], stdout=f).wait()
        f.close()

        if euc_proc != 0:
            sys.stderr.write(f"WARNING: Euclidean compact script returned non-zero exit code for {pdbfile} {euc_proc}!\n")
        else:
            print(f"Processed: {pdbfile}")


def get_secondary_structure():
    """
    Runs DSSP to generate secondary structure data for input files.
    Results for each model are saved in a new directory: 'OUTPUT_DIR/dssp'
    """
    dssp_outpath = assert_directory(OUTPUT_DIR + "/dssp")
    mkdssp = join(QA_DIR,"tools/mkdssp")
    dssp2dataset = join(QA_DIR,"scripts/dssp2dataset.pl")
    dssp2secondary = join(QA_DIR,"scripts/dssp2secondary.pl")
    cmd = "perl %s %s %s %s %s" % (dssp2secondary, PDB_DIR, mkdssp, dssp2dataset, dssp_outpath)
    cmd_args = cmd.split()

    print(cmd)
    dssp_proc = subprocess.Popen(cmd_args).wait()

    if dssp_proc != 0:
        sys.stderr.write(f"WARNING: There was a problem running DSSP (code: {dssp_proc}), please check the output files at {dssp_outpath}\n")

    return dssp_outpath


def gen_surface(dssp_outpath):
    """
    Generates a Surface Score for each input file.
    Results are saved as surface_score.txt
    """
    surface_outfile = join(OUTPUT_DIR, "surface_score.txt")
    cmd = "perl %s %s %s" % (QA_DIR+"/scripts/surface_score.pl", dssp_outpath, surface_outfile)
    cmd_args = cmd.split()

    print(cmd)
    surface_proc = subprocess.Popen(cmd_args).wait()

    if surface_proc != 0:
        sys.stderr.write(f"WARNING: There was a problem generating surface score (code: {surface_proc}), please check the output at {surface_outfile}\n")


def gen_wea(dssp_outpath):
    """
    Calculates Weighted Exposed Area.
    Results are saved as weighted_exposed_area.txt
    """
    outfile = join(OUTPUT_DIR, "weighted_exposed_area.txt")
    cmd = "perl %s %s %s" % (QA_DIR+"/scripts/weight_exposed.pl", dssp_outpath, outfile)
    cmd_args = cmd.split()

    print(cmd)
    wea_proc = subprocess.Popen(cmd_args).wait()

    if wea_proc != 0:
        sys.stderr.write(f"WARNING: There was a problem generating weighted exposed area (code: {wea_proc}), please check the output at {outfile}\n")


def gen_sbrod():
    """
    Generates scores for the following using SBROD: Residues, Hydrogen Bonds, Solvation, and Backbone.
    Results are saved under 'OUTPUT_DIR/sbrod' (4 files):
        backboneatom.sbrod
        hbonds.sbrod
        residues.sbrod
        solvation.sbrod
    """
    output_path = join(OUTPUT_DIR, "sbrod")
    bins_folder = join(OUTPUT_DIR, "sbrod_bin")
    pdb_folder_name = PDB_DIR.split('/')[-1]
    script_path = join(QA_DIR, "tools/SBROD/protein_scoring/scripts/")

    features_cmd = "./generate_dataset.sh %d %s %s" % (1, PDB_DIR, bins_folder)
    move_cmd = "mv %s/binaries/%s %s" % (bins_folder, pdb_folder_name, output_path)

    cwd = os.getcwd()
    os.chdir(script_path)   # change working directory
    os.system(features_cmd) # run generate_dataset.sh
    if len(os.listdir(f"{bins_folder}/binaries/{pdb_folder_name}")) != 4:
        sys.stderr.write(f"WARNING! Some SBROD features might have failed to generate! ({pdb_folder_name})\n")
    os.chdir(cwd)           # change back to original directory
    os.system(move_cmd)     # move generated files to output folder
    assert not os.path.isdir(f"{bins_folder}/binaries/{pdb_folder_name}") # should be empty after moving files
    os.system("rm -r %s" % bins_folder) # sbrod generates a folder we don't need, gotta delete it


def gen_topQA():
    """
    Uses TopQA to generate an 8x8x8 representation of the protein's topology, reshaped
    to a 1-D array (512 columns total).
    Results are saved at: 'OUTPUT_DIR/topqa/TopQAFeatures.txt'
    """
    topqa_output = OUTPUT_DIR + "/topqa"
    topqa_script = QA_DIR + "/tools/TopQA/script/TopQA_pred.py"
    pdb2img_path = QA_DIR + "/tools/TopQA/resources/pdb2img"
    model_path = QA_DIR + "/tools/TopQA/resources/Model_520.01"
    cmd = "python %s %s %s %s %s" % (topqa_script, pdb2img_path, model_path, PDB_DIR, topqa_output)
    cmd_args = cmd.split()

    print(cmd)
    topqa_proc = subprocess.Popen(cmd_args).wait()

    if topqa_proc != 0:
        sys.stderr.write(f"WARNING: There was a problem generating TopQA features (code: {topqa_proc}), please check the output at {topqa_output}\n")


def gen_ANDIS():
    """
    Runs ANDIS and reads the scores it generates.
    Final output is saved as andis.txt in the output directory.
    """
    cwd = os.getcwd() # save current working directory before we change it
    andis_folder = ANDIS_PATH

    # Create a pdblist as required by ANDIS
    pdblist_path = PDB_DIR + "/pdblist.txt"
    with open(pdblist_path, 'wt') as pdblist:
        for pdb in os.listdir(PDB_DIR):
            if not "pdblist" in pdb:
                pdblist.write(pdb+"\n")

    # Run ANDIS and direct output to a temporary file (tmp_file)
    tmp_file = OUTPUT_DIR + "/andis.tmp"
    andis_cmd = "./ANDIS %s %s" % (PDB_DIR, pdblist_path.split('/')[-1]) # don't use full path for pdblist
    os.chdir(andis_folder) # change working directory

    print(andis_cmd)
    f = open(tmp_file, 'w')
    andis_proc = subprocess.Popen(andis_cmd.split(), stdout=f).wait() # run ANDIS
    f.close()

    if andis_proc != 0:
        sys.stderr.write(f"\nWARNING: ANDIS returned non-zero exit code ({andis_proc}), please check the output file ({tmp_file}) and your ANDIS installation.\n")

    os.chdir(cwd)          # return to original working directory
    os.remove(pdblist_path) # delete the pdblist file

    # Read data from tmp file, a dictionary is used to avoid accidental duplicates
    data = {}
    with open(tmp_file, 'rt') as tmp:
        for line in tmp:
            if '*' in line: continue # skip non-data lines
            if ('warning' in line.lower()) or ('error' in line.lower()):
                sys.stderr.write(f"ERROR: There was a problem running ANDIS! Please check the output file ({tmp_file}) and your ANDIS installation.\n")
                sys.exit()
            name, value = line.split()[0].split('.')[0], line.split()[1]
            data[name] = value
    if andis_proc == 0: os.remove(tmp_file)

    # Write results to file
    with open(OUTPUT_DIR+"/andis.txt", 'wt') as outfile:
        for name in data:
            line = name + ' ' + data[name] + '\n'
            outfile.write(line)


def gen_energy_terms():
    """
    PyRosetta Energies (this will generate 12 of them)
    Results are saved as model_name.rosetta under 'OUTPUT_DIR/rosetta'.
    """
    outpath = assert_directory(OUTPUT_DIR+"/rosetta")
    log_folder = assert_directory(OUTPUT_DIR + "/log")
    f = open(log_folder+"/rosetta.log", 'w') # output log file
    ros_proc = subprocess.Popen(["python3", QA_DIR+"/scripts/ros_energy.py",
                                 "-d", PDB_DIR,
                                 "-o", outpath], stdout=f).wait()
    f.close()

    if ros_proc != 0:
        sys.stderr.write("WARNING: There was a problem generating rosetta energies (Code: %d)! Please check the output and PyRosetta installation.\n" % ros_proc)




################################################################################
# READ OUTPUTS -----------------------------------------------------------------
################################################################################

def read_score_file(filepath, feature_data, label):
    """
    Reads from an output (e.g. surface_score.txt) file that contains data for each PDB.
        filepath: path to the score file
        feature_data: the feature data dictionary (see generate_features_file())
        label: label of feature in the feature_data dictionary to append data to

    Returns the modified feature_data dictionary containing the new data
    """
    if not os.path.isfile(filepath):
        sys.stderr.write(f"ERROR: {filepath} is missing. Please check your logs, input files, installation, etc. and try running {sys.argv[0]} again.\n")
        sys.exit()

    with open(filepath, 'rt') as file:
        for line in file:
            name, score = line.split()[0].split('.')[0], line.split()[1]
            feature_data[name][label] = score
    return feature_data


def read_topqa_output(topqa_pred_file, features, feature_labels):
    """
    Reads from the TopQA output file and appends the data to variables used by generate_features_file().
    Each TopQA feature label starts with "TQA" followed by a number corresponding to its indice (TQA1, TQA2, ...)
        topqa_pred_file: path to the TopQAFeatures.txt
        features: The feature data dictionary from generate_features_file()
        feature_labels: The list of feature labels from generate_features_file()
    Returns the modified features and feature_labels containing TopQA data
    """
    with open(topqa_pred_file, 'rt') as tqfile:
        for line in tqfile:
            cols = line.split() # get a list of the row's columns
            name = cols[0].split('.')[0] # exclude the file extension when grabbing the name
            for i in range(len(cols[1:])): # for each topqa feature
                col_name = 'TQA'+str(i+1)
                if not (col_name in feature_labels):
                    feature_labels += [col_name] # add to feature labels list if they haven't already
                assert col_name in feature_labels # just to make sure
                features[name][col_name] = cols[i+1]
    return features, feature_labels


def read_rosetta_file(ros_file):
    """
    Reads an rosetta output file and returns the data in the form of a list containing
    sublists of data.

    For example: [ [row1_column1, row2_column1],
                   [row1_column2, row2_column2],
                   [row1_column3, row2_column3] ]
    """
    ros_data = [[0] for _ in range(12)] # instantiate a list of empty lists [ [] , [], ..., [] ]
    with open(ros_file, 'rt') as rfile:
        for line in rfile:
            columns = line.split()
            for i in range(12): # for each column
                ros_data[i].append( float(columns[i+1]) ) # append each column to the corresponding list in ros_data
    return ros_data


def generate_features_file():
    """
    Reads all the generated files and puts the data together into a single, combined file called
    "feature_data.txt"

    All data in this output file is raw (not normalized). Normalized data will
    be output to a seperate file by another function.

    Missing features are automatically handled at the end of the function. A pdb that has
    a value missing for a feature will have that value set to 'inf'. PDBs with 'inf' values
    are skipped in the normalization and prediction step.
    """
    feature_labels = [ # codenames for each of the features used (topqa and rosetta energies are appended to this list down below)
                "EUCL",      # euclidean compact score
                "SURF",      # surface score
                "WEXA",      # weighted exposed area
                "RESIDUES",  # sbrod feature 1
                "HBONDS",    # sbrod feature 2
                "SOLVATION", # sbrod feature 3
                "BACKBONE",  # sbrod feature 4
                "ANDIS"     # ANDIS output
                ]
    energy_terms = [
                "R_ENV",    # residue environment
                "R_PAIR",   # residue pair interactions
                "R_CBETA",  # c_beta density
                "R_VDW",    # steric repulsion
                "R_RG",     # radius of gyration
                "R_CENPACK",# packing
                "R_CO",     # contact order
                "R_HS",     # statistical potentials
                "R_SS",
                "R_RSIGMA",
                "R_SHEET",
                "R_CEN_HB"  # centroid hydrogen bonding
    ]
    feature_labels += energy_terms

    features = {} # this will store all the data we will write to a file

    # Instantiate empty dictionaries for each pdb file
    pdbFilenames = os.listdir(PDB_DIR)
    for pdb in pdbFilenames:
        features[pdb.split('.')[0]] = {}

    # 1. Read score files (results for all the PDBs are in one file)
    ### Read Surface Score ###
    ##########################
    features = read_score_file(OUTPUT_DIR+"/surface_score.txt", features, "SURF")

    ### Read Weighted Exposed Area ###
    ##################################
    features = read_score_file(OUTPUT_DIR+"/weighted_exposed_area.txt", features, "WEXA")

    ### SBROD Features ###
    ######################
    features = read_score_file(OUTPUT_DIR+"/sbrod/residues.sbrod", features, "RESIDUES")
    features = read_score_file(OUTPUT_DIR+"/sbrod/hbonds.sbrod", features, "HBONDS")
    features = read_score_file(OUTPUT_DIR+"/sbrod/solvation.sbrod", features, "SOLVATION")
    features = read_score_file(OUTPUT_DIR+"/sbrod/backboneatom.sbrod", features, "BACKBONE")

    ### ANDIS ###
    #############
    features = read_score_file(OUTPUT_DIR+"/andis.txt", features, "ANDIS")

    ### TopQA Features ###
    ######################
    #topqa_pred_file = join(OUTPUT_DIR, "topqa/TopQAFeatures.txt")
    #assert os.path.isfile(topqa_pred_file), "Couldn't find %s!" % topqa_pred_file
    #features, feature_labels = read_topqa_output(topqa_pred_file, features, feature_labels)

    # 2. Read features that each PDB has their own file for
    for pdb in pdbFilenames:
        name = pdb.split('.')[0] # name of the model without the file extension

        ### Read Euclidean Compact Scores ###
        #####################################
        euc_file = join(OUTPUT_DIR, "euclidean", name+".euc")
        if not os.path.isfile(euc_file):
            sys.stderr.write(f"WARNING! Couldn't find {euc_file}!\n")
            # TODO: Placeholder value? Might be handled below already
        else:
            f = open(euc_file, 'rt')
            euc_score = next(f).split()[1]
            features[name]["EUCL"] = euc_score
            f.close()

        ### ROSETTA Energy Scores (Averages) ###
        ########################################
        # env pair cbeta vdw rg cenpack co hs ss rsigma sheet cen_hb
        ros_file = join(OUTPUT_DIR, "rosetta", name+".rosetta")
        if not os.path.isfile(ros_file):
            print(f"WARNING! Couldn't find {ros_file}!")
            # TODO: See todo for Euclidean Compact above
        else:
            ros_data = read_rosetta_file(ros_file)
            # Calculate the averages
            i = 0
            for term_label in energy_terms:
                features[name][term_label] = str( sum(ros_data[i])/len(ros_data[i]) )
                i += 1

        ### WRITE FEATURES TO FILE ###
        ##############################
        # TODO: .replace(' ', '_')
        col_width = 25 # columns are padded with spaces so that they line up, this is the max width.
        feat_dir = assert_directory(OUTPUT_DIR + "/features")
        features_filepath = join(feat_dir, name + ".feat")
        with open(features_filepath, 'w+') as feature_file:
            header = "#NAME" + " "*(col_width-5)
            row = name + " "*(col_width-len(name)) # instantiate row data with model name
            for label in feature_labels:
                header += label + " "*(col_width-len(label))
                if label in features[name]:
                    row += features[name][label] + " "*(col_width-len(features[name][label]))
                else: # if the feature is missing due to an error, etc.
                    if label == "ANDIS": # TODO: this is a temporary fix for ANDIS failing on large proteins; please find a better imputation method
                        row += "0" + " "*(col_width-1)
                    else:
                        row += str( float('inf') ) + " "*(col_width-3)
            feature_file.write(header + '\n')
            feature_file.write(row + '\n')

    ### COMBINE FEATURE FILES TOGETHER ###
    ######################################
    feat_dir = join(OUTPUT_DIR, "features")
    with open(OUTPUT_DIR + "/feature_data.txt", 'wt') as combined_feature_file:
        # Write the header
        header = "#NAME" + " "*(col_width-5)
        for label in feature_labels:
            header += label + " "*(col_width-len(label))
        combined_feature_file.write(header + '\n')

        # Write each row of data, skipping rows with feature values missing ('inf')
        print('\n')
        for filename in os.listdir(feat_dir):
            feat_file = open(feat_dir + '/' + filename, 'rt')
            next(feat_file) # skips the header
            row = next(feat_file)

            if str( float('inf') ) in row:
                print(f"SKIPPING: {filename} (some features failed to generate).")
            else:
                combined_feature_file.write(row)
            feat_file.close()


def normalize_data(feature_filepath, output_file):
    """
    Takes feature_data.txt and normalizes each of the features appropriately. It will
    also filter out some rosetta features that will not be used.
    This function will output the results as "feature_data_normalized.txt"
    """
    with open(feature_filepath, 'rt') as raw_data_file:
        header = next(raw_data_file)[1:] # grab the header, excluding the '#' at the beginning
        #first_topqa_column = header.split().index("TQA1")
        Model = NamedTuple("Model", header.split())
        print(Model._fields, '\nNumber of Columns: ', len(Model._fields))

        # Write the new header
        new_header = "#NAME TARGET GDT EUCL SURF WEXA RESIDUES HBONDS SOLVATION BACKBONE ANDIS R_ENV R_PAIR R_CBETA R_VDW R_CENPACK R_CEN_HB "
        #for i in range(len(header.split()[first_topqa_column:])):
            #new_header += "TQA{} ".format(i+1)
        new_header += '\n'

        new_data = ""
        feature_vector = []
        topqa_vector = []
        pdb_names = []

        # Append each row of normalized data to the new_data variable
        for row in raw_data_file:
            row_data = row.split()
            row_data = [row_data[0]] + [float(i) for i in row_data[1:]] # first column is a string, convert the rest to floats
            m = Model._make(row_data) # convert data to a named tuple
            assert len(row_data) == len(Model._fields), "Could not parse the following:\n%s" % str(row_data) # make sure the columns match up

            # Put normalized features together into a new list
            feat = [m.EUCL, m.SURF, m.WEXA,
             rescale(m.RESIDUES, -52.879900, 149.769823),
             rescale(m.HBONDS, -100.174316, 392.383050),
             rescale(m.SOLVATION, -713.281763, 36156.863968),
             rescale(m.BACKBONE, -254.847080, 840.384436),
             rescale(m.ANDIS, -19255.470000, 474.210000),
             sigmoid(m.R_ENV),
             rescale(m.R_PAIR, -0.485115, 0.458073),
             rescale(m.R_CBETA, 0.330233, 0.557829),
             sigmoid(m.R_VDW),
             rescale(m.R_CENPACK, -0.332316, 1.444480),
             sigmoid(m.R_CEN_HB)]

            assert len(feat) == 14

            # Concatenate the features together into a single string and append to new_data
            row_string = m.NAME + " N/A N/A " # gdt and target name are used for training, but have no use for scoring
            for f in feat:
                row_string += str(f) + ' '
            #for f in row_data[first_topqa_column:]:
                #row_string += str(f) + ' '
            row_string += '\n'
            new_data += row_string

            # Append features to vectors / lists
            feature_vector.append(feat)
            #topqa_vector.append( row_data[first_topqa_column:] )
            pdb_names.append(m.NAME)

    with open(output_file, 'w+') as outfile:
        outfile.write(new_header)
        outfile.write(new_data)
    print("Normalized features have been saved as: ", output_file)
    return np.array(feature_vector), np.array(topqa_vector), pdb_names # return these variables so we can directly use them for scoring







################################################################################
# FEATURE PAIRS LAYER ----------------------------------------------------------
################################################################################
def concat_features(feature_set1, feature_set2):
    """
    Appends feature_set1 to the end of feature_set2.
    Both variables are a list of arrays.

    e.g. feature_set1 = [ [1,2,3], [4,5,6] ]
         feature_set2 = [ [0,0,0], [1,1,1] ]
         concat_features(feature_set1, feature_set2) returns [ [0,0,0,1,2,3], [1,1,1,4,5,6] ]
    """
    final_array = []
    for layer1_features, layer2_features in zip(feature_set1, feature_set2):
        combined_feat = np.concatenate( (layer2_features, layer1_features) )
        final_array.append(combined_feat)
    return np.array(final_array)


def gen_layer2_features(feature_data, models):
    """ Makes a prediction using every possible feature pair for one row of data. """
    features = [] # should be 91 features
    i = 0
    for index1 in range(NUM_FEATURES):
        for index2 in range(index1+1, NUM_FEATURES):
            pred = models[i].predict([[ feature_data[index1], feature_data[index2] ]])[0]
            features.append(pred)
            i += 1
    assert len(features) == 91
    return features


def get_layer2_features(data, models):
  """
  Generates the second feature layer that contains predictions made by
  all the trained classifiers on every row of data.

  data: a list containing nested lists with feature vectors for each pdb (each list represents a row of data).
    e.g. [ [0.512, 0.194, ..., 0.914], [0.102, ...] ]

  i.e. this funciton calls gen_layer2_features for each row of data and appends the result to a cumulative list.
  """
  layer2_features = []
  for feature_inputs in data:
    layer2_features.append(gen_layer2_features(feature_inputs, models))
  return layer2_features







###################################
###################################




def main():
    global PDB_DIR

    print("\n###############################" +
          "\n# GENERATING FEATURES... ------" +
          "\n###############################")
    assert_directory(OUTPUT_DIR)
    tmp_folder = assert_directory(OUTPUT_DIR+"/tmp")
    os.system("cp -r %s %s" % (PDB_DIR, tmp_folder)) # copy files to output folder to resolve permission issues
    PDB_DIR = tmp_folder + '/' + PDB_DIR.split('/')[-1]


    print("\nGenerating Secondary Structure...")
    dssp_outpath = get_secondary_structure()
    print("Secondary Structure Finished!\n")

    print("\nGenerating Euclidean Compact Score...")
    gen_euclidean()
    print("Euclidean Finished!\n")

    print("\nCalculating Surface Score...")
    gen_surface(dssp_outpath)
    print("Surface Score Finished!\n")

    print("\nCalculating Weighted Exposed Area...")
    gen_wea(dssp_outpath)
    print("Weighted Exposed Area Finished!\n")

    print("\nCalculating Rosetta Energies...")
    gen_energy_terms()
    print("Rosetta Finished!")

    print("\nGenerating SBROD features...")
    gen_sbrod()
    print("SBROD Finished!\n")

    print("\nRunning ANDIS...")
    gen_ANDIS()
    print("ANDIS Finished!\n")

    #print("\nGenerating TopQA features...")
    #gen_topQA()
    #print("TopQA Finished!\n")

    print("\nGenerating features file...")
    generate_features_file()
    print("Done!\n")


    # Create a seperate file containing normalized data
    layer1_features, topqa_features, pdb_names = normalize_data(join(OUTPUT_DIR, "feature_data.txt"), join(OUTPUT_DIR, "feature_data_normalized.txt"))

    ### Load Pickled Models ###
    ###########################
    pair_models = load(join(QA_DIR, "models", "SVR_Pairs_1.model"))
    top1_clf = load(join(QA_DIR, "models", "SVR_Top1.model"))
    ave_clf = load(join(QA_DIR, "models", "SVR_Averages.model"))
    print("\nPairs (layer 1):", pair_models[0])
    print("Final (TOP1):", top1_clf)
    print("Final (AVE):", ave_clf)

    # Generate layer 2 features (2-feature models / feature pairs)
    layer2_features = get_layer2_features(layer1_features, pair_models)
    layer2_features = concat_features(layer1_features, layer2_features) # append the original features
    #layer2_features = concat_features(topqa_features, layer2_features) # append the TopQA features
    print(layer2_features, layer2_features.shape)

    print("\n###############################" +
          "\n# SCORING MODELS... -----------" +
          "\n###############################")
    try:
        scores = ave_clf.predict(layer2_features)
        top1_scores = top1_clf.predict(layer2_features)
    except ValueError: # should occur if there's only one sample (one input pdb)
        layer2_features = layer2_features.reshape(1, -1)
        scores = ave_clf.predict(layer2_features)
        top1_scores = top1_clf.predict(layer2_features)

    # Make sure they're all between 0 and 1 (most of them should be)
    for i in range(scores.shape[0]):
        scores[i] = 0.0001 if scores[i] < 0 else scores[i]
        scores[i] = 0.9999 if scores[i] > 1 else scores[i]
        top1_scores[i] = 0.0001 if top1_scores[i] < 0 else top1_scores[i]
        top1_scores[i] = 0.9999 if top1_scores[i] > 1 else top1_scores[i]
    #print(scores, '\n', pdb_names)

    # Find the best model predicted by the "top 1" model and move it to the top of the "average" model's rankings
    top1_scores = [(name, score) for name, score in zip(pdb_names, top1_scores)]
    top1_scores = sorted(top1_scores, key = lambda x: x[1])[::-1] # from greatest to smallest
    best_model = top1_scores[0][0]
    #print("\nTop1 Classifier:\n", top1_scores, '\n')

    # Update score in "average" model predictions to be the max score plus 0.00001
    for i in range(scores.shape[0]):
        if pdb_names[i] == best_model: # name matches
            scores[i] = scores.max() + 0.00001
            print(f"\nNOTICE: Set {best_model} as the highest-scoring model.\n")
            break

    # Sort results
    results = [(name, score) for name, score in zip(pdb_names, scores)]
    results = sorted(results, key = lambda x: x[1])[::-1] # from greatest to smallest
    #print(results)

    # Write to file
    predictions_filename = "_predictions.txt"
    with open(join(OUTPUT_DIR, predictions_filename), 'w+') as pred_file:
        for r in results:
            row = f"{r[0]} {r[1]}\n"
            print(row)
            pred_file.write(row)

    # Clear the tmp folder
    os.system("rm -r %s" % OUTPUT_DIR+"/tmp/"+PDB_DIR.split('/')[-1])

    # Finish
    end_time = time.monotonic()
    total_time = timedelta(seconds=end_time - start_time)
    print("\n--------------------------------------------------------------")
    print(f"Finished! Predictions have been saved as {predictions_filename}")
    print("Duration: %s" % str(total_time))
    print(f"Date: {time.asctime()}")
    print("User: %s" % pwd.getpwuid( os.getuid() )[ 0 ])



if "__main__" in __name__:
    main()
