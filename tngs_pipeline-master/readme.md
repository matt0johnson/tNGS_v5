# **Snakemake targeted NGS data analysis pipeline**

This analysis pipeline is designed to align fastq files to a reference genome and produce SNV/indel/CNV variant calls over targeted intervals using the Snakemake workflow management system. The steps undertaken are described in the snakefiles subfolder, and are chained together as follows:

![alt tag](https://github.com/rdemolgen/tngs_pipeline/blob/master/workflow.png?raw=true)

## Usage

### **1. Activate virtual environment:**

        source /mnt/Data4/targeted_sequencing/virtual_environ/bin/activate

### **2. Clone pipeline Github repository:**

        git clone https://rdemolgen@github.com/rdemolgen/tngs_pipeline.git <workbatch_folder>    

### **3. Update sample_list_all.csv file in workbatch/snakefiles folder with information from lab excel sheet:** 

### **4. Edit snakefiles/config.yaml:**

	Enter workbatch, ftp, flowcell, base_path, batch_id

### **5. Run snakemake command:**

        nohup snakemake -k --cores 16 2>&1 &

### **6. Deactivate virtual environment:**

        deactivate

