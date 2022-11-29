# Stability Analysis Tool Documentation

## Installation
------------

1.  Download [Docker Desktop](https://www.docker.com/products/docker-desktop/)
2.  Set up [WSL](https://learn.microsoft.com/en-us/windows/wsl/setup/environment)
    *   Restart your machine after the install
    *   Open Ubuntu using the start menu and create a Linux username & password
3.  Download [this](https://drive.google.com/drive/u/0/folders/1o5gMu7jSA6pBcJ_qEH3ZPyKTo_ute4iy) folder to somewhere on your local drive where you want the software tool to live
4.  Install the Docker image
    1.  Open Command Prompt
    2.  Navigate to the Stability\_Tool folder you just downloaded
        *   Use `pwd` to see your current directory, `ls` to see the contents of the directory, and `cd [directory]` to switch to a given directory
    3.  Run `docker pull pymesh/pymesh`
        *   You should see pymesh/pymesh appear in the Images section in the Docker Desktop app
    4.  Run `docker build -t pymesh/visualization .`
        *   You should see pymesh/visualization appear in the Images section in the Docker Desktop app

  

## Running the Tool
----------------

1.  To use the stability analysis tool, you first need to enter the Docker image you created:
    *   Run `docker run --name stability_tool_container -it --rm -v "$(pwd):/root" pymesh/visualization bash`
        *   This will link your current directory to `/root` in the Docker container
2.  Run the tool: `python3 main.py [arguments]`
3.  Arguments:
    *   `-f [filename]` or `--filename [filename]` : filename of input STL (mandatory)
    *   `-m [mass]` or `--mass [mass]` : mass of the buoy in kilograms (mandatory)
    *   `-c [X Y Z]` or `--cg [X Y Z]` : location of buoy center of gravity in millimeters (mandatory)
    *   `-r [resolution]` or `--resolution [resolution]` : angular resolution of analysis in degrees (optional, defaults to 45 degrees)
    *   `-a [accuracy]` or `--accuracy [accuracy]` : Buoyancy accuracy of water line computation (in kilograms). (Optional, default is 1 gram accuracy. i.e. buoyancy will match mass within 1 grams.)
    *   `-v [axis]` or `--rotation [axis]` : rotation axis to spin the buoy about (X or Y) (optional, defaults to X)
4.  Example: run `python3 main.py -f buoy.stl -m 0.5 -c 0 0 50 -r 5 -a 0.001 -v y` to analyze the file "buoy.stl", asserting a mass of 0.5 kg, a CG location of (0, 0, 50) mm, an angular resolution of 5 degrees, a mass accuracy of 0.001 kg, and revolving the mesh about the Y axis
5.  When done with the tool, run `exit` to leave the Docker container

  

## Notes
-----

*   The mesh file must be saved in the `Stability_Tool/inputs` directory
*   The output of the analysis will be saved to a folder in the `Stability_Tool/outputs` directory
    *   The output folder will contain a plot of the righting moments (positive is self-righting), a CSV containing relevant metrics, and an animation of the buoy
*   The mesh must have units of millimeters
*   Frequently, Solidworks creates [non-manifold](https://www.instructables.com/Non-manifolds-Your-Worst-3D-Printing-Nightmare/) meshes when exporting STLs. If the stability tool fails (`Input mesh is not PWN!`), you should fix the mesh with a mesh fixing software. [This one](https://www.formware.co/onlinestlrepair) has worked ok.
*   The default "Up" axis in Solidworks is Y, while the up axis in the stability tool is Z. You'll likely have to rotate your part before using the tool.
