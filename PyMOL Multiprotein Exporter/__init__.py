#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: guymayneord
"""

from __future__ import absolute_import
from __future__ import print_function
from pymol import cmd
from pymol.Qt import QtWidgets, QtCore
from pymol.Qt.utils import loadUi
from os import path
from platform import system as current_system
from subprocess import Popen, call
from time import sleep


# =============================================================================
# Parameters
# =============================================================================
SURFACE_TYPES = ["Surface", "Cartoon",
                 "Sticks", "Mesh",
                 "Lines", "Ribbon",
                 "Wire", "Licorice",
                 "Dots", "Spheres"]

# After X seconds close the dialog.
COUNTDOWN = 3

Dialog = None


# =============================================================================
# PyMOL start-up
# =============================================================================
# Standard function for adding pymol plugin
def __init_plugin__(app=None):
    from pymol.plugins import addmenuitemqt
    addmenuitemqt("Multiprotein Exporter", runPluginGUI)


def runPluginGUI():
    global Dialog
    if Dialog is None:
        Main_app = MainWindow()
        Dialog = Main_app.returnDialog()
    Dialog.show()


# =============================================================================
# Exporting function:
# =============================================================================
class ExportingThread(QtCore.QThread):
    prog_update = QtCore.pyqtSignal(dict)

    def __init__(
            self, dir_out, export_type, output_obj_list, prog_bar_obj,
            progress_label):

        QtCore.QThread.__init__(self, None)
        self.dir_out = dir_out
        self.export_type = export_type
        self.output_obj_list = output_obj_list
        self.progress_bar_obj = prog_bar_obj
        self.progress_label = progress_label

    def run(self):
        # Hide everything in the scene so it can be individually displayed.
        cmd.hide('everything')
        print('==============================================')
        print('Creating surfaces and exporting...:')
        print('==============================================')

        # Show each object as determined by the user.
        for obj_no in range(len(self.output_obj_list)):
            cmd.disable('all')
            obj_name = self.output_obj_list[obj_no].text()
            print(f"    -{obj_name}")
            # Send over the data to the progress bar to reflect updates
            self.prog_update.emit({"Job_Number": obj_no + 1,
                                   "Total_jobs": len(self.output_obj_list),
                                   "Object_name": obj_name})

            cmd.show_as(self.export_type.lower(), obj_name)
            cmd.enable(obj_name)
            output_name = path.join(self.dir_out, f"{obj_name}.wrl")
            cmd.save(output_name, ['wrl'])
            self.output_obj_list[obj_no].setStyleSheet(
                "background-color: green")

        print(' ')
        self.prog_update.emit({"Job_Number": len(self.output_obj_list),
                               "Total_jobs": len(self.output_obj_list),
                               "Object_name": None,
                               "Last_file": output_name})

        print('==============================================')
        cmd.disable('all')
        print('Complete')
        print('==============================================')


# =============================================================================
# GUI:
# =============================================================================
class MainWindow(QtWidgets.QMainWindow):
    def __init__(self, *args, **kwargs):
        super(MainWindow, self).__init__(*args, **kwargs)
        uifile = path.join(path.dirname(__file__), 'Main_menu.ui')
        self.Gui = loadUi(uifile, Dialog)
        self.scroll_widget_area = QtWidgets.QWidget()
        self.v_layout = QtWidgets.QVBoxLayout()
        self.scroll_widget_area.setLayout(self.v_layout)
        self.Gui.scrollArea.setWidget(self.scroll_widget_area)

        self.individual_obj_checkboxes = []
        self.connectInterfaceButtons()
        self.populateObjAction()

    def connectInterfaceButtons(self):
        self.Gui.refresh_button.clicked.connect(
            self.populateObjAction)

        self.Gui.select_all.clicked.connect(
            lambda: self.selectAllAction(select_opt=True))

        self.Gui.select_none.clicked.connect(
            lambda: self.selectAllAction(select_opt=False))

        self.Gui.browse_button.clicked.connect(
            lambda: self.browseButtonAction(self.Gui.output_directory))

        self.Gui.Submit_button.clicked.connect(self.submitButtonAction)
        self.Gui.Cancel_button.clicked.connect(self.closeButtonAction)
        self.populateExportsAction(self.Gui.export_type)

    def returnDialog(self):
        # Just returns the dialog for handling in PyMOL
        return self.Gui

    def populateExportsAction(self, dropdown_box_to_update):
        for each_surface in SURFACE_TYPES:
            dropdown_box_to_update.addItem(each_surface)

    def selectAllAction(self, select_opt):
        for each_checkbox in self.individual_obj_checkboxes:
            each_checkbox.setChecked(select_opt)

    def populateObjAction(self):
        # Destroy previous list but maintain reference
        for i in range(len(self.individual_obj_checkboxes)):
            self.individual_obj_checkboxes[0].deleteLater()
            del self.individual_obj_checkboxes[0]

        # Make new list of new objects.
        pyMOL_objects_list = cmd.get_names('objects', 0, '(all)')
        for each_object in pyMOL_objects_list:
            new_checkbox = QtWidgets.QCheckBox(each_object)
            new_checkbox.setChecked(True)
            self.v_layout.addWidget(new_checkbox)
            self.individual_obj_checkboxes.append(new_checkbox)

    def submitButtonAction(self):
        if self.Gui.output_directory.text() != "":
            obj_for_export = []
            for each_checkbox in self.individual_obj_checkboxes:
                if each_checkbox.isChecked():
                    obj_for_export.append(each_checkbox)

            if len(obj_for_export) == 0:
                print('==============================================')
                print('No objects selected in scene.')
                print('==============================================')
            else:
                print(' ')
                print('==============================================')
                print(f"{len(obj_for_export)} objects selected, exporting all to {self.Gui.output_directory.text()}")
                print('==============================================')
            self.Exporting_thread = ExportingThread(
                self.Gui.output_directory.text(),
                self.Gui.export_type.currentText(),
                obj_for_export,
                self.Gui.progress_bar,
                self.Gui.progess_text)

            self.Exporting_thread.prog_update.connect(self.updateProgress)
            self.Exporting_thread.start()
        else:
            self.Gui.output_directory.setStyleSheet("background: red;")

    def updateProgress(self, recieved_dict):
        progress_percent = int((recieved_dict["Job_Number"] / recieved_dict["Total_jobs"]) * 100)
        self.Gui.progress_bar.setValue(progress_percent)

        if recieved_dict["Object_name"] is not None:
            self.Gui.progess_text.setText(f"Exporting job {recieved_dict['Job_Number']} of {recieved_dict['Total_jobs']}: {recieved_dict['Object_name']}")
        else:
            self.Gui.progess_text.setText("Complete!")
            # Open the directory we've just written the files to
            self.openFileInExplorer(recieved_dict["Last_file"])

            for i in range(COUNTDOWN, 0, -1):
                self.Gui.progess_text.setText(f"Closing dialog in {i-1}")
                QtWidgets.QApplication.processEvents()
                sleep(1)
            self.closeButtonAction()

    def browseButtonAction(self, textbox_to_update):
        textbox_to_update.setStyleSheet("background: None;")
        # Get the directory location
        option = QtWidgets.QFileDialog.ShowDirsOnly
        option |= QtWidgets.QFileDialog.DontUseNativeDialog
        filename = str(QtWidgets.QFileDialog.getExistingDirectory(self,
                                                                  "Select Directory for output",
                                                                  options=option))

        if filename:
            textbox_to_update.setText(filename)

    def closeButtonAction(self):
        # Need to clear up the interface once a job is complete.
        self.Gui.output_directory.setText("")
        self.Gui.progress_bar.setValue(0)
        self.Gui.progess_text.setText("")
        for i in range(len(self.individual_obj_checkboxes)):
            self.individual_obj_checkboxes[0].deleteLater()
            del self.individual_obj_checkboxes[0]
        self.Gui.close()

    def openFileInExplorer(self, file_path):
        if current_system() == 'Windows':
            Popen(r'explorer /select,"' + f'{file_path}"')
        if current_system() == 'Darwin':
            call(['open', '-R', '-aFinder', f"{file_path}/"])
