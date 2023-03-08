# Ubuntu 16.04
To compile C code:   
make  
To compile Python wrapper  
swig -python -c++ -extranative s6fits.i  
make _s6fits.so  

# Ubuntu 20.04
## Requirements
* Please Make sure the following libraries are installed:  
    ```
    sudo apt install swig  
    sudo apt install gnuplot
    sudo apt install imagemagic
    ```
* [Miniconda](https://docs.conda.io/en/latest/miniconda.html) is recommended here.
## Compile the s6fits library
1. create and activate conda environment
    ```
        conda create -n seti python=3.9
        conda activate seti
    ```
1. modify Makefile.py39
    ```
    CFLAGS=-I$(IDIR) -I/home/$USER/.conda/envs/seti/include/python3.9 -Wall -fPIC -g
    LIBS= -L /home/$USER/.conda/envs/seti/lib
    ```
    **Note:** If you use conda, and your conda env name is "seti", you don't need to change anything.
2. convert c methods to python methods
    ```
    swig -c++ -python s6fits.i
    make -f Makefile.py39 _s6fits.so
    ```
3. copy the three files to your data directory
    ```
    cp _s6fits.so /data/serendip6_data
    cp s6fits.py /data/serendip6_data
    cp s6_plot_waterfall_py39.py/serendip6_data
    ```
## Plot waterfall
* You should be able to plot the waterfall figure now!
    ```
    cd /data/serendip6_data
    python s6_plot_waterfall_py39.py xxxx.fits
    ```
    **Note:** When you run the script, please make sure you're using the same python environment when you compile the s6fits library.
* you can use `display` to look at the figure.
    ```
    disply xxx.wf.png
    ```
    **NOTE:** Please use `ssh -XY xx@xxx` to log into your server.
## Error and fix
*   When you are trying to plot the waterfall figure, you may have the issue like this:
    ```
    ImportError: /lib/x86_64-linux-gnu/libp11-kit.so.0: undefined symbol: ffi_type_pointer, version LIBFFI_BASE_7.0
    ```
* If so, please cd to your python lib directory, and rename libffi.so.7.   
Then you need to create a new symbolic:
    ```
    cd /home/$USER/.conda/envs/seti/lib
    mv libffi.so.7 libffi.so.7.bak
    ln -s /lib/x86_64-linux-gnu/libffi.so.7 libffi.so.7
    ```
