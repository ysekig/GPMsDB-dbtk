# GPMsDB-dbtk

GPMsDB-dbtk v1.0.1 was released on March 7, 2023. 

GPMsDB-tk/GPMsDB-dbtk are software toolkits for assigning taxonomic identification to user-provided MALDI-TOF mass spectrometry profiles obtained from bacterial and archaeal cultured isolates. They take advantages of a newly developed database of protein mass profiles predicted from ~200,000 bacterial and archaeal genome sequences. This toolkit is also designed to work with customized databases, allowing microbial identification based on user-provided genome/metagenome-assembled genome (MAG) sequences. The GPMsDB-dbtk is open source and released under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International Public License. 

GPMsDB-dbtk is used for customizing the GPMsDB with user-provided genomes and MAGs/SAGs. 

Please post questions and issues related to GPMsDB-dbtk on the Issues section of the GitHub repository.

## Installing and using GPMsDB-dbtk

Prerequisites
* Python (version 3.7 or higher)
* Cython (version 0.29.1 or higher)
* [GPMsDB-tk](https://github.com/ysekig/GPMsDB-tk) (version 1.0.1 or higher) and relevant database ([R01-RS95](https://zenodo.org/record/7703483#.ZAbPNS_3Jf0))
* biolib (https://github.com/dparks1134/biolib) (version 0.1.8 or higher)
* matplotlib (version 3.5.0 or higher)
* Prodigal (https://github.com/hyattpd/Prodigal) (version 2.6.3 or higher)
* HMMER (https://github.com/EddyRivasLab/hmmer) (version 3.3.2 or higher)

In the source directory, the following command will compile and install the software in your python environment.
```bash
git clone https://github.com/ysekig/GPMsDB-dbtk
cd GPMsDB-dbtk
python setup.py install
```
During the installation, you may see some deprecation warnings like “easy_install command is deprecated” but this will not cause any issues for GPMsDB-dbtk.

GPMsDB-dbtk requires an environment variable named GPMsDB_PATH to be set to the directory containing the unarchived reference data.
```bash
export GPMsDB_PATH=/path/to/release/package/
```

### Features

* Genome(s) to massDB:
  * genome_wf     -> Full genomes to ms data workflow
  * list_db       -> List genome entries in the custom ms database
  * update_db     -> Add peak_list(s) to the custom ms database
  * remove_genome -> Delete entries from the custom ms database
			
## Bug Reports

Please report bugs through the GitHub issues system, or contact Yuji Sekiguchi (y.sekiguchi@aist.go.jp)

## Copyright

Copyright (C) 2023 Yuji Sekiguchi, National Institute of Advanced Industrial Science and Technology (AIST)

This package (the majority of the scripts in the package) is under the conditions of the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International Public License. See LICENSE for further details.
