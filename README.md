# PyMOL-Multiprotein-Exporter
PyMOL plugin to export each object in a PyMOL scene to a individual .wrl file. 

The purpose of this add-on was to allow low-resource computers to export objects in PyMOL to .wrl format, a common 3D object file. 
For standard single-protein sessions, this may not have much application, however for large multi-protein sessions, this will aid in exporting different objects in the scene. 

The project consists of a directory containing a .py and .ui file, and a zip file which is the compressed version of the directory. The add-on can be installed either by importing the .py file (the .ui file is automatically recocognised in the installation), or the .zip file. If there are any issues importing via the .py file, feel free to create your own .zip file from the directory for safety.

When the add-on is opened, the user can select which objects in the scene are to be exported, and the type of export (Cartoon, surface etc.) and the directory to output the files to. 
Once the export is complete, the dialog will close and the newly exported files will be displayed.

