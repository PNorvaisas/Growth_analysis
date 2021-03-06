# Growth analysis
Scripts for Biolog (`Biolog.py`) and other bacterial growth (`Justgrowth.py`) experiment analysis in multiwell plates.
- Calculates bacterial growth parameters such as integral (AUC) and growth rate
- Plots curves for growth, growth in log scale and growth rate (dOD/dt)
- Exports *tidy* data as time series, summary and experiement description
## Configuration
### MacOS
- Install [Xcode](https://developer.apple.com/xcode/) from App Store and open it to accept license (might ask to update your OS to most recent version).
- Install [Homebrew](https://brew.sh) by running in Terminal:
   ```bash
   /usr/bin/ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"
   ```
- Install git version control system using brew:
    ```bash
    brew install git
    ```
- Install Python using brew:
   ```bash
   brew install python
   ```
- Configure your Python:
   ```bash
   brew unlink python && brew link python;
   brew doctor
   ```
   > In some cases might require mannually setting default Python to Homebrew version:
   > ```bash
   > echo "export PATH=/usr/local/bin:/usr/bin:$PATH" >> ~/.bash_profile;
   > cp /usr/local/bin/python2 /usr/local/bin/python
   >```
- Install necessary Python packages:
   ```bash
   pip install numpy, scipy, matplotlib
   ```
- Clone `Growth_analysis` GitHub repository into selected `<location>`:
   ```bash
   cd <location>;
   git clone git@github.com:PNorvaisas/Growth_analysis.git
   ```
- Make sure scripts are executable by running:
   ```bash
   chmod +x *.py
   ```
- Add `<location>` to your PATH:
   ```bash
   echo "export PATH=<location>:$PATH" >> .bash_profile
   ```
- Log Out & Log In
- Test scripts by running in any directory:
    ```bash
    Biolog.py -h
    Justgrowth.py -h
    ```
- Script should ask you to install the necessary packages that are missing.
- If you can't run these scripts in any directory, there might be a problem with the way you have set PATH. To test this, try running scripts in the `<location>` where you have cloned git repository. If you still can't do that, there might be other problems.

### Windows
- Install [Anaconda](https://conda.io/docs/user-guide/install/windows.html)
    - Choose Anaconda instead of Miniconda (reduced functionality)
    - You need to add Anaconda to your PATH!
- Install [git](https://git-scm.com/download/win) for Windows
- [Navigate](https://www.computerhope.com/issues/chusedos.htm) to a selected `<location>` using `Command Prompt`
- Clone `Growth_analysis` GitHub repository into a selected `<location>`:
   ```batch
   cd <location>;
   git clone git@github.com:PNorvaisas/Growth_analysis.git
   ```
- In this `<location>` try running:
    ```batch
    python Biolog.py -h
    python Justgrowth.py -h
    ```
- If everything is correct, scripts should automatically check for missing packages and install them.
---
## Biolog.py
- `Biolog.py` uses `Design.xlsx` and `Biolog_metabolites.csv` tables to generate output to a selected folder (defaults to `Output`). 
- `Design` (template: `Biolog_Design_template.xlsx`) file can be any table (`.xls`, `.xlsx`, `.csv`) that contains columns:
    - **File** - name of the plate reader file.
    - **Plate** - name of the Biolog plate as present in the `Biolog_metabolites_EcoCycU.csv` table.
    - **Reader** - name of the reader used to collect data: **Tecan** or **Biotek**.
    - **Strain** - bacterial strain used in an experiment.
    - **Type** - experiment type: **Control**/**Treatment**. Defines **line color** in generated plots.
    - **Replicate** - biological replicate of the experiment. Defines **linetype** in generated plots.
    - Any additional variables of interest which can be applied to the whole plate, like **Replicate**, **Drug**, etc.
- As a rule of a thumb, experiment description in `Design` file should be as detailed as possible.
### Running in MacOS, Unix:
- Navigate to the working directory (where your plate reader and `Design.xlsx` files are) in `Terminal` and run:  
    ```bash
    Biolog.py -i Design.xlsx -d Biolog_metabolites_EcoCycU.csv -o Output_dir
    ```
- To get information about input arguments:
    ```batch
    Biolog.py -h
    ```
### Running in Windows
- Copy both `Biolog.py` and `Biolog_metabolites_EcoCycU.csv` to where your plate reader and `Design.xlsx` files are.
- [Navigate](https://www.computerhope.com/issues/chusedos.htm) to the working directory (where your plate reader and `Design.xlsx` files are) in `Command Prompt` and run: 
    ```batch
    python Biolog.py -i Design.xlsx -d Biolog_metabolites_EcoCycU.csv -o Output_dir
    ```
- To get information about input arguments:
    ```batch
    python Biolog.py -h
    ```
---
## Justgrowth.py
- `Justgrowth.py` uses `Design.xlsx` and `Pattern.xlsx`  tables to generate output to a selected folder (defaults to `Output`). 
- `Design` (template: `Justgrowth_Design_template.xlsx`) file can be any table (`.xls`, `.xlsx`, `.csv`) that contains columns:
    - **File** - names of the plate reader files.
    - **Pattern** - names of `Pattern` tables, which contain information about the variables investigated in a particular plate
    - **Reader** - name of the reader used to collect data: **Tecan** or **Biotek**
    - Any additional variables of interest which can be applied to the whole plate, like **Replicate**, **Date**, etc.
- `Pattern` (template `Justgrowth_Pattern_template.xlsx`) file can be any **Excel** table (`.xls`, `.xlsx`) with multiple sheets, each of which describes different variable in the plate and particular value assigned to each well:
    - **Strain**, **Drug**, **Drug_conc**, etc.
- As a rule of a thumb, experiment description in `Design` file should be as detailed as possible
- `Design` file can be used to describe variables which are **shared** in the whole plate.
- `Pattern` file can be used to describe variable which are **unique** to each well in the plate.
### Runnig in MacOS, Unix
- Navigate to the working directory (where your plate reader and `Design.xlsx` files are) in `Terminal` and run: 
    ```batch
    Justgrowth.py -i Design.xlsx -o Output_dir
    ```
- To get information about input arguments:
    ```batch
    Justrowth.py -h
    ```
### Windows
- Copy `Justgrowth.py` to where your plate reader and `Design.xlsx` files are.
- [Navigate](https://www.computerhope.com/issues/chusedos.htm) to the working directory (where your plate reader and `Design.xlsx` files are) in `Command Prompt` and run: 
    ```batch
    python Justgrowth.py -i Design.xlsx -o Output_dir
    ```
- To get information about input arguments:
    ```batch
    python Justgrowth.py -h
    ```
---
## Growth.py
- `Growth.py` is a script that combines the functionality of both `Justgrowth.py` and `Biolog.py`. It **uses the same input files** as the previous scripts and **generates the same output**. Given that most of the functions are shared between `Justgrowth.py` and `Biolog.py`, it only makes sense to join them into one, more consistent framework.
- `Growth.py` uses `growth` or `biolog` flags at the end of the command argument to choose `Justgrowth.py` or `Biolog.py` functionality correspondingly.
- `Growth.py` also has `full` flag, which generates more extensive output in terms of data tables and figures.
### Runnig in MacOS, Unix
- Navigate to the working directory (where your plate reader and `Design.xlsx` files are) in `Terminal`: 
-  For **regular reads** which are used with `Pattern` files (like `Justgrowth.py`):
    ```batch
    Growth.py -i Design.xlsx -o Output_dir growth
    ```
    For more extensive output:
    ```batch
    Growth.py -i Design.xlsx -o Output_dir growth full
    ```
-  For **Biolog** reads which are used with `Biolog_meabolites.csv` table (like `Biolog.py`):
    ```batch
    Growth.py -i Design.xlsx -m Biolog_metabolites.csv -o Output_dir biolog
    ```
    For more extensive output:
    ```batch
    Growth.py -i Design.xlsx -m Biolog_metabolites.csv -o Output_dir <biolog/growth> full
    ```
- To get information about input arguments:
    ```batch
    Growth.py -h
    ```
### Windows
- Copy `Growth.py` to where your plate reader and `Design.xlsx` files are.
- [Navigate](https://www.computerhope.com/issues/chusedos.htm) to the working directory (where your plate reader and `Design.xlsx` files are) in `Command Prompt`: 
-  For **regular reads** which are used with `Pattern` files (like `Justgrowth.py`):
    ```batch
    python Growth.py -i Design.xlsx -o Output_dir growth
    ```
    For more extensive output:
    ```batch
    python Growth.py -i Design.xlsx -o Output_dir growth full
    ```
-  For **Biolog** reads which are used with `Biolog_meabolites.csv` table (like `Biolog.py`):
    ```batch
    python Growth.py -i Design.xlsx -m Biolog_metabolites.csv -o Output_dir biolog
    ```
    For more extensive output:
    ```batch
    python Growth.py -i Design.xlsx -m Biolog_metabolites.csv -o Output_dir <biolog/growth> full
    ```
- To get information about input arguments:
    ```batch
    python Growth.py -h
    ```
---
# Output
- Description to be added..

