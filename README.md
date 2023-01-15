# PyMOL-Multiprotein-Exporter
PyMOL plugin to export each object in a PyMOL scene to a individual .wrl file.

The purpose of this add-on was to allow low-resource computers to export objects in PyMOL to .wrl format, a common 3D object file. 
For standard single-protein sessions, this may not have much application, however for large multi-protein sessions, this will aid in exporting different objects in the scene. 

When the add-on is opened, the user can select which objects in the scene are to be exported, and the type of export (Cartoon, surface etc.) and the directory to output the files to. 
Once the export is complete, the dialog will close and the newly exported files will be displayed.
