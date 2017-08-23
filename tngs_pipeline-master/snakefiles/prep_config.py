import subprocess, csv, os, fnmatch, yaml

# Open config.yaml file
with open('config.yaml', 'r') as f:
    config  = yaml.load(f)

# Get required parameters from config file
ftp = config["ftp"]
flowcell = config["flowcell"]
workbatch = config["workbatch"]
path = workbatch + "/rawdata/"

# Download fastq and md5sum files
command = 'wget -r -nd -nv -A fastq.gz %s/%s/Sample_v* -P %s --ignore-case' % (ftp, flowcell, path)
command += '; wget -r -nd -nv -A fastq.gz.md5sum %s/%s/Sample_v* -P %s --ignore-case' % (ftp, flowcell, path)
command += '; cd %s ; md5sum -c *.md5sum > md5sum.out 2>&1' % (path)
command += ';chmod 444 *'
subprocess.run(command, shell=True)

# Set empty samples array

samples = []

# Collect lists of fastq files

r1 = []
r2 = []

for file in os.listdir(path):
    if not file.startswith('.'):
        if fnmatch.fnmatch(file, "*R1*.fastq.gz"):
            r1.append(file)
        elif fnmatch.fnmatch(file, "*R2*.fastq.gz"):
            r2.append(file)

# Define sample class

class Sample:
    def __init__(self, sample_name, gender, library, panel):
        self.sample_name = sample_name
        self.gender = gender.upper()
        self.library = library
        self.panel = panel

    def add_forward(self, sample_name):
        for fastq in r1:
            if sample_name in str(fastq):
                self.forward = fastq

    def add_reverse(self, sample_name):
        for fastq in r2:
            if sample_name in str(fastq):
                self.reverse = fastq

# Create samples from csv input file

with open("sample_list_all.csv") as infile:
    csvreader = csv.reader(infile, delimiter = ',')
    # Skip header row
    row2 = next(infile).rstrip()
    # For all sample rows, add sample data to Sample class
    for row in csvreader:
        row[3] = Sample(row[1], row[4], str.split(row[5])[0], str.split(row[5])[1])
        row[3].add_forward(row[1])
        row[3].add_reverse(row[1])
        samples.append(row[3])

# Output to config.yaml

with open("config.yaml", "a") as outfile:    
    for sample in samples:
        outfile.write('  %s : [%s,%s,%s,%s,%s]\n' % (sample.sample_name, sample.forward, sample.reverse, sample.gender, sample.library, sample.panel))
