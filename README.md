# NMAT: Nano-micromotor analysis tool

Tool to perform the most typical statistical analysis of nano- and micromotors, including MSD, MSAD, trajectory plotting, autocorrelation velocity, and performs different types of fittings and motion parameter extraction (speed, diffusion coefficient, rotational diffusion time, rotational speed), including averages between particles of the same conditions.

## Citation and contact

If you use this script for your research, please cite it as follows, adapting if necessary (check lisence file):

> The analysis of motion was performed with the Python-based Nano-micromotor Analysis Tool (NMAT) v. *X* (https://github.com/rafamestre/NMAT-nanomicromotor-analysis-tool). 

where *X* is your version. Feel free to add the GitHub address as a reference, footnote or between brackets.

If you want your research to be hightlighted here as a user of this tool, please email me at r.mestre[at]soton[dot]ac[dot]uk. For issues, errors or requests for changes/collaborations, email me. 

### Tracking of motion

If you don't know how to track your particles' motion, I recommend to use my Nano-micromotor Tracking Tool ([NMTT](https://github.com/rafamestre/NMTT-nanomicromotor-tracking-tool)), which can track particles or other objects from videos of bright-field or fluorescent microscopy. With this tool, you can track the position of one or more of your particles at the same time using machine learning techniques, and the results include tracking videos and position *vs* time results that can be used directly with the NMAT to analyse the motion. 

## Installing dependencies

It's recommended to use Python 3.6 or higher, although this script is backwards compatible to Python 2. It has been tested in all three Windows, Linux and MacOS.


In any case, very few not built-in modules are necessary (numpy, matplotlib, seaborn, scipy) and the script checks if they exist before launching for the first time, and installs them otherwise. So technically it shouldn't even be necessary to install them from the requirements.txt file. 

## How to run

#### Option 1: Python already installed

If you have Python already installed in your computer (not with an Anaconda distribution, or apart from it), you should be able to run the file by double-clicking in it. Double click on the file should open a graphic user interface (GUI) to start the analysis. If double click doesn't work, check that the file opens with Python by right-clicking in it and selecting "open with" and then "Python". 

If that doesn't work, you can try to run it from the terminal by writing:

```
python NMATv1.py
```

(Assuming the NMATv1.py file is in the home directory, otherwise write [pathtofile]/NMATv1.py).

#### Option 2: with Anaconda through terminal

If the previous doesn't work, it's better to use [Anaconda](https://www.anaconda.com/products/individual) to run this code and install the necessary packages. Once Anaconda is installed, open an Anaconda prompt and type:

```
python NMATv1.py
```

(Assuming the NMATv1.py file is in the home directory, otherwise write [pathtofile]/NMATv1.py). The first time you run it, some packages not built-in with Python are installed for you to use them. It might take some seconds and some stuff might be written in the prompt, but eventually the GUI will appear. Next time you open it, those packages won't be installed and it'll be much faster.

#### Option 3: with Anaconda with a different environment

If that didn't work neither, you might have some packages in your base environment that clash with the script. In this case it's better to create a separate conda environment. 

After installing Anaconda, open your Anaconda prompt and write:

```
conda create --name NMAT python=3.6
```

This will create an environment with the Python version 3.6, which works for sure. Once the environment has been created, activate it by typing:

```
conda activate NMAT
```

You should be able to then run the script by typing:

```
python NMATv1.py
```

(Assuming the NMATv1.py file is in the home directory, otherwise write [pathtofile]/NMATv1.py). As before, the first time run it, some packages might be installed and it might take some time.

#### Option 4: with Anaconda through Spyder

You might choose to run it through an Integrated Development Environment (IED) like Spyder. To install it, write in the Anaconda prompt (in whichever environment you choose):

```
conda install spyder
```

Once it's been installed, you can open it by simply writing in the prompt:

```
spyder
```

Once opened, you can drag the *.py* file into the window or open it in File\Open. Clicking F5 or the button "Run file" will launch the script.

#### Manually installing packages

If the automatic installation of packages gives any problem, you can install them by yourself by downloading the requirements.txt file in your home directory and typing:

```
pip install -r requirements.txt
```

(Note: it is assumed the requirements.txt file is in the home directory, otherwise you should write "[pathtofile]/requirements.txt"



## Instructions to run

A Graphic User Interface (GUI) like the following one appears:


<p float="left">
  <img src="https://user-images.githubusercontent.com/13152269/148577649-335d2db1-222e-42cc-8136-824907ffa0d8.png" width="250" />
</p>

### Opening files

By clicking in **Browse** a pop-up window to select a folder will appear. Select the folder where you data is. **Important note**: it is recommended to create separate folder for different conditions and not mix data that you wouldn't aggregate because they represent different experimental conditions. For instance, you might want to create 5 folders for 5 different fuel concentrations. Inside each folder, **only** add results from particles that used that fuel concentration. This is because the script can perform averages of particles, and you wouldn't want to mix different conditions and perform an average.

### Data distribution

Each file needs to have data about the X position, Y position and time in seconds. The script accepts data that's organised in two different ways:

#### Horizontal

You might choose to have your files with the data horizontally. In this way, you can have one single file that has information about many particles, as in the following screeshot (from the [NMTT](https://github.com/rafamestre/NMTT-nanomicromotor-tracking-tool) tool):

![image](https://user-images.githubusercontent.com/13152269/148578627-5547e375-bdd8-4abe-86dc-2ae3be235dc4.png)

In this example, each particle is introduced by the word "Particle" with the particle number on the right. The script will read whatever's below "Particle" and look for the time, X and Y data. It's case insensitive and as long as "time", "x" and "y" appear, it should be able to detect the rows. It doens't matter if you write *time (seconds)*, *time (s)* or simply *time* - as long as the word "time" appears, it should be fine. Likewise for X and Y, you can write *x* or *X* or *x (micrometers)* or *x position*, as long as there's an "x", it will detect the row of data. It is quite flexible, but also not invincible, so try to stick to simple labels and avoid other type of rows to be present in the file. For instance, if there's one row that's named "Velocity", since there's a "y" in it, the script might confuse it by the Y position.

X and Y are assumed to be in micrometers and don't need to start at 0 um. Time is assumed to be in seconds and doesn't need to start at 0 s.

Introduce each new particle with a new "Particle" instance. The horizontal distribution of data is especially useful if you don't want each particle to be in a separate file. 

#### Vertical

If each one of your tracked particles is in a different file, you can organise it vertically, as in the following screenshot:

![image](https://user-images.githubusercontent.com/13152269/148580103-8e086176-54f7-47c1-aa1c-714a192ac9a1.png)

It's worth stressing that for data ordered vertically **each file should contain data for a single particle**. The first row can be "Particle" **or not**. If "particle" appears in the first row, the script will assign that particle label to it. If no particle label was given, it will assign the label 0.

The next (or first) row is the header. As before, the script will look for "time", "x" and "y", and will not assume any position. Again, you can write *time (s)* or simply *time* or *TiMe* and the script will understand it. After this, it will read the data vertically, until the file is finished.

### File type and delimeter

The NMAT accepts two types of files, .csv or .txt files. Likewise, it will also accept two delimiteres, either tab (\t) or commas (,). Obviously, for .csv files, you should select the .csv type and the , delimeter. For .txt files, you can choose. 

### File reading

Once you have browsed to your desired folder and have selected the type of file, delimeter and data distribution, the NMAT will automatically look for files. It will do this dynamically, so if your data is distributed vertically and you switch from "Vertical" to "Horizontal" version and viceversa, the script will see and lose the files correspondingly. 

The NMAT tells you how many files were found of the specified type and inside of those how many valid particles it found. A valid particle is one that was introduced by a particle label (except for the vertical type) and had correct "Time", "X" and "Y" data. If you used horizontal distribution of data, you might have more particles than files. For instance, if you had 10 files with 5 particles each, you'll have 10 files and 50 particles. But if you used vertical distribution, since each file can contain only information about a single particle, you can only have as many particles as files in the folder.

**Important**: do not leave more .csv or .txt files in the folder than you need. The script is fairly robust and should not read particles if they don't contain all the information (time, x, y), but there might be ways to get around that. 

### Analysis

First of all, all the analyses can be made **individually** (by default), but also **averaging between particles** (highly encouraged). Select the option "Do averages between particles?" below to do this. In each folder, the script will output the averaged results but also will create a folder called "Individual" with results for each single particle.

Below the reading options, there are the analysis options:

Analysis             |  Description
:-------------------------:|:-------------------------:
MSD | Select to calculate the Mean Square Displacement (MSD) of each particle
Trajectory | Select to plot the trajectory of each particle in micrometers. The script creates individual plots of trajectories that **have the same scale**, therefore all the particles are comparable. You can use a vector graphic program like Inkscape to open the .svg files and place the trajectories together to compare
Instantaneous velocity | Select to plot the instantaneous velocity of the particles. **Note**: use only this for very ballistic trajectories, and even then use it with care. Instantaneous velocity is not a good measure of the actual speed of the nanomotors, especially if they are showing only enhanced diffusion. Check our [paper](https://arxiv.org/abs/2007.15316) for more info
MSAD | Select to calculate the Mean Square Angular Displacement (MSAD) of each particle, which can be useful to find active rotation, e.g., particles with non-zero rotational speed (chiral). Chiral particles have a quadratic MSAD. Check our [paper](https://arxiv.org/abs/2007.15316) for more info
Velocity autocorrelation | Select to calculate the velocity autocorrelation of each particle. This is useful to find ballistic behaviour or directionality in the particles. If the autocorrelation decreases in an exponential manner with time, there might be directionality in the motion (a quadratic MSD should be visible). If it goes immediately to zero except at *t = 0*, the motion is mainly diffusive

This is the basic analysis, mainly for visualisation purposes, but there are some options to extract even more parameters from this:

#### Local exponent of MSD

Select to calculate the local exponent of the MSD. This is the derivative of the logarithm of the MSD with respect to the logarithm of time, or the exponent ![formula](https://render.githubusercontent.com/render/math?math=\alpha) assuming an MSD of the form ![formula](https://render.githubusercontent.com/render/math?math=MSD(t)=t^\alpha). This is useful to see if the motion is purely diffusive (![formula](https://render.githubusercontent.com/render/math?math=\alpha=1)), purely ballistic (![formula](https://render.githubusercontent.com/render/math?math=\alpha=2)) or is between regimes (![formula](https://render.githubusercontent.com/render/math?math=1<\alpha<2)

#### Fit MSD

Select to perform different types of fittings to the MSD. You can choose between Linear, Quadratic, Third Order, Fourth Order, Full Equation or All. The election of type of fitting should come from the type of motion that is observed and the physical assumptions behind it. 

Are your particles very small and clearly show enhanced diffusion? Use a linear fitting and only report enhanced diffusion, never speed (not even instantaneous).

Are your particles clearly showing directional motion (or you're not sure but they could)? Use a quadratic fitting, or even better, use third or fourth order fittings, which greatly reduce the error associated to the parameters. In this case **do not report the diffusion value**: this parameter is greatly overfitted and can have extreme errors. If you use third or fourth order fittings, you might also obtain an approximation for the rotational diffusion time, although it might have some error. Second order fitting works better when the rotational diffusion time is very far away from the fitting interval. If the rotational diffusion time is inside or close to your fitting interval, third or fourth order fittings help improve the fitting, since a quadratic equation is not the best approximation in this scenario.

Are you very sure your particles are showing classical active motion, without deviations? Then you might want to fit it directly to the Full Equation, as presented in the classical paper by [Howse *et al.* (2007)](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.99.048102). This equation is complex and therefore the fitting might fail very easily or be unstable, so unless you have a large number of particles, the results might not be so good. It's recommended to do polynomial fittings as explained above.

Are you unsure of everything? Just click All and do all the fittings (or plot only the MSD and see how it looks like). However, be aware of [HARKing](https://en.wikipedia.org/wiki/HARKing) or hypohtesising after your results are known, or [p-hacking](https://en.wikipedia.org/wiki/Data_dredging), choosing the type of fitting or parameters that give you the best results. Your choice of fitting should come from your *a priori* assumptions of the physical phenomena taking place.

A very detailed and thorough discussion of the types of fittings and the errors associated to them can be found in our [paper](https://arxiv.org/abs/2007.15316).

#### Fit MSAD

Select to perform different fittings to the MSAD. You can only select Linear or Quadratic fitting, or Both. A particle shows a linear MSAD if their motion is non-chiral, where the slope if the rotational diffusion coefficient. If there is chiral motion, the MSAD will look quadratic, and by performing a quadratic fitting you can extract the rotational speed of the particle.

#### Logarighmic fitting of the MSD

Apart from the typical polynomial fittings, you can perform a linear fitting to the logarithm of the MSD, that gives you information about the local exponent, assuming a form ![formula](https://render.githubusercontent.com/render/math?math=MSD(t)=t^\alpha). It's similar to the calculation of the local exponent, but in this case it's not local, it's global. This works better for shorter fitting intervals where the shape of the MSD hasn't changed much. It's useful when the particle is showing strange motion and the assumptions of the classical active Brownian motion approximations might not hold (e.g. you cannot assign a speed that's accompanying the quadratic term in the MSD). In this case, you might just want to compare the exponent of the MSD: a larger exponent means more active motion, whatever type it might be. 


### Interval of interest

By default, when you have at least 1 valid particle in your folder, the interval of interest is automatically set to 1/10th of the total duration of the largest trajectory. This is a convention in the field, since the MSD is known to behave unreliably (with high standard deviation) at long times. You can however change this interval if you wish to. Typical intervals of interest are between 1 and 4 s, but this largely depends on the trajectory length and the number of particles. This is the interval that will be used to perform the fittings explained above.

Most of the analyses are shown twice, one with a "longer" version, and another with a "shorter" version. For instance, for transparency's sake, the MSD is shown with its long version, which is very unreliable at long times and might show oscillations, but also with its short version, which is cut at the end of the interval of interest. This is the version that should be used to plot the results or perform any fittings. Remember that the shorter the interval, the more reliable the results, unless you have very low FPS and very few points in the interval.

### Only trajectories longer than...

It doesn't matter if you have particles of approximately 30 s of trajectory, but then you have some other particles with only 5 s. You should not discard that data, because it's still useful. In most experimental cases, your system of nano- or micromotors will be ergodic, meaning that averaging in time or making an ensemble average are identical. That is, it would be the same having 10 particles with 10 s trajectories separately or concatenating all those trajectories to form a single particle with a trajectory of 100 s. So if you have a particle of 30 s and a short one of 5 s, it's the same as having a single particle of 35 s. This is only true if you're **averaging between particles and in time**. Obviously, if you're considering independent trajectories, the 30 s one will be more realiable than the 5 s one because it's longer and you have more data points. In any case, it's better to try and keep the same time and make it always 10 times larger than the interval of interest (at least 20 s if possible). More on this in our [paper](https://arxiv.org/abs/2007.15316).

However, you should eliminate particles that fall below your interval of interest. So if you want to look at particles from 0 s to 5 s and you have a trajectory of 4 s, you should elimiate it. This is because when going from 4 s to 4.1 s, there'll be a jump in the MSD because of this discontinuity, and will make it look weird. Therefore, in "Only trajectories longer than ..." you should write at least the same number as the interval of interest, although feel free to restrict it even more and write even a longer time if you don't really trust your short videos. The information above with the "number of valid particles" will update itself when you select/deselect this option and when you write the numbers (make sure it's not showing "no valid particles" because you wrote 0 s!).

## Results

As mentioned above, each type of analysis will create a folder. Inside that folder, if you selected "Do averages" (highly encouraged), the average result will be there. Then, inside the "Individual" folder, the results for each particle will also be there, sometimes in long and short version. There are three types of files:

File             |  Description
:-------------------------:|:-------------------------:
Plot in .png | The plot of the analysis in .png format. If it's a fitting, it will contain the parameter estimation in the legend. 
Plot in .svg | The plot of the analysis in .svg format. This is vector graphic format, which can be modified with a program like [Inkscape](https://inkscape.org/), Adobe Illustrator or similar. With this programs you can directly change the colors, widths, etc., without losing quality and use the plots directly in your paper.
Data in .csv | The data in .csv format, if you want to plot yourself the results. For individual particles, the first row is always the name of the file that contained that specific particle, for instance *\user\tracking1.csv*. All particles are named subsequently, starting for 0, independently of the particle label. In order to differentiate them, in the second row, you will have the particle label. If the file corresponds to the particle label 23 inside the file *\user\tracking1.csv* (that contained many particles), you can find it easily here. Maybe the results of that particle looked really strange and you want to check out the video closely to figure out what's going on!
Summary file | Some folders will have a summary file (mainly fittings). It will contain all the individiual fitting results with the average value, standard deviation and standard error of the mean.

### Note on the fittings:

You can chose to extract the parameters from the fittings in two ways: 

1) doing individual fittings and then making the average of the parameters;

2) doing the average of the MSD and then one single fitting. 

Both techniques should, in most scenearios, yield approximately the same results, given the ergodicity of the system. However, there are some things to take into account:

* If some of the trajectories are very short, you should not use technique 1). What I mentioned above in subsection "Only trajectories longer than..." is valid when making ensemble averages (between all particles). In this case, you would be extracting parameters from very short trajectories, which would be unreliable, and then making an average with the rest of fittings. In this case, it's better to use technique 2), because you'll be using the information about the short trajectories to get a better estimation of the results. In this case, the error of the parameters it's not a simple standard deviation from making an average, but **the error of the fitting**, which is reported by the NMAT in the file. This error, however, doesn't tell you how anything about the distribution of, let's say, speed in your particles, only about how good was your fitting (if you had enough particles, it'll be a very good fitting for sure). If your particles have very different velocities and you're interested in this distribution, this is not the correct approach. If you're only interested in the overall effect, it should be fine. You also won't be able to perform statistical tests on your data to compare between conditions.
* If your trajectories are stable and long, you can choose either technique 1) or 2), it should really make no difference. Technique 1) will give you an estimate of the parameters with the standard error of the mean, while technique 2) will give you the error of the fitting. Actually, if you want to see how "stable" your results are, go and compare the individiual fitting summary file with the average fitting and see if they give roughly the same value. 



