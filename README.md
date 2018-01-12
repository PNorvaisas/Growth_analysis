# Growth analysis
Scripts for bacterial growth and Biolog experiment analysis in multiwell plates.
- Calculates bacterial growth parameters such as integral (AUC) and growth rate
- Plots curves for growth, growth in log scale and growth rate (dOD/dt)
- Exports *tidy-ish* data as time series and summary
## Configure
### MacOS
- Install [Xcode](https://developer.apple.com/xcode/) from App Store and open it to accept license (might ask to update your OS to most recent version).
- Install [Homebrew](https://brew.sh) by running in Terminal:
   ```bash
   /usr/bin/ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"
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
- Clone [Growth_analysis](git@github.com:PNorvaisas/Growth_analysis.git) GitHub repository into selected `<location>`:
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
- Restart

    