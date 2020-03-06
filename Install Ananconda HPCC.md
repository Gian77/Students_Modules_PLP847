# to run Python 2 and Biopython in the HPCC cluster, Michigan State University 
# you need the anaconda platform. Let's get started on how to get your own anaconda!

*Very useful links*

https://wiki.hpcc.msu.edu/display/~longnany@msu.edu/Using+conda<br/>
http://docs.anaconda.com/anaconda/install/linux/

Download Anaconda here https://www.anaconda.com/download/#linux

Upload the installer in the HPCC and run

`wget https://repo.anaconda.com/archive/Anaconda2-5.3.0-Linux-x86_64.sh`<br/>
`bash ~/Anaconda2-5.3.0-Linux-x86_64.sh`<br/>

Review the license agreement, clic Yes to agree.
Enter Yes to repend the Anaconda2 install location to PATH in your /home/<user>/.bashrc
The installer finishes and displays Thank you for installing Anaconda2
To make the installation effective run

`source ~/.bashrc`

After your install is complete, verify it by opening Anaconda Navigator (this is optional)

`anaconda-navigator`

Check the insalled version by running

`conda --version`

To use Biopython we need to create an environment with python 2 and then 
install biopython 

`conda create --name py2 python 2.7`

You will probably get prompted to install new packages, type ‘y’ and hit "Enter". 
When it is done, display your list of environments with

`conda info --envs` <br/>

Look in your newly created environment if you have the packages you wanted

`conda activate py2` <br/>
`conda list` <br/>

Now you can add a package to the created environment

`conda install -n py2 biopython`

To use conda and activate one of your environments and use it on the terminal
or inside a job type

`conda activate py2`<br/>
`<your commands>`<br/>
`conda deactivate`

# You are done. Enjoy your Biopython on HPCC @ MSU!

If you want to add another tool/software to one of your environemnts,
for eample, PEAR https://cme.h-its.org/exelixis/web/software/pear/
an ultrafast, memory-efficient and highly accurate pair-end read merger,
you can do the following:

`conda install -c bioconda pear`

Now. if you want to use PEAR, you just need to load the
py2 conda environment.
