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
surface_types = ["Surface", "Cartoon",
                 "Sticks", "Mesh",
                 "Lines", "Ribbon",
                 "Wire", "Licorice",
                 "Dots", "Spheres"]

# =============================================================================
# PyMOL start-up
# =============================================================================
def __init_plugin__(app=None):
    from pymol.plugins import addmenuitemqt
    addmenuitemqt("Multisubunit exporter", run_plugin_gui)


dialog = None


def run_plugin_gui():
    global dialog
    if dialog is None:
        main_app = MainWindow()
        dialog = main_app.return_dialog()
    dialog.show()

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
            self.prog_update.emit({"Job_Number": obj_no+1,
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
# GUI element:
# =============================================================================
class MainWindow(QtWidgets.QMainWindow):
    def __init__(self, *args, **kwargs):
        super(MainWindow, self).__init__(*args, **kwargs)
        uifile = path.join(path.dirname(__file__), 'Main_menu.ui')
        self.ui_obj = loadUi(uifile, dialog)
        self.scroll_widget_area = QtWidgets.QWidget()
        self.v_layout = QtWidgets.QVBoxLayout()
        self.scroll_widget_area.setLayout(self.v_layout)
        self.ui_obj.scrollArea.setWidget(self.scroll_widget_area)

        self.individual_obj_checkboxes = []
        self.connect_interface_buttons()
        self.populate_obj_action()

    def connect_interface_buttons(self):
        self.ui_obj.refresh_button.clicked.connect(
            self.populate_obj_action)

        self.ui_obj.select_all.clicked.connect(
            lambda: self.select_all_action(select_opt=True))

        self.ui_obj.select_none.clicked.connect(
            lambda: self.select_all_action(select_opt=False))

        self.ui_obj.browse_button.clicked.connect(
            lambda: self.browse_button_action(self.ui_obj.output_directory))

        self.ui_obj.Submit_button.clicked.connect(self.submit_button_action)
        self.ui_obj.Cancel_button.clicked.connect(self.close_button_action)
        self.populate_exports_action(self.ui_obj.export_type)

    def return_dialog(self):
        # Just returns the dialog for handling in PyMOL
        return self.ui_obj

    def populate_exports_action(self, dropdown_box_to_update):
        for each_surface in surface_types:
            dropdown_box_to_update.addItem(each_surface)

    def select_all_action(self, select_opt):
        for each_checkbox in self.individual_obj_checkboxes:
            each_checkbox.setChecked(select_opt)

    def populate_obj_action(self):
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

    def submit_button_action(self):
        if self.ui_obj.output_directory.text() != "":
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
                print(f"{len(obj_for_export)} objects selected, exporting all to {self.ui_obj.output_directory.text()}")
                print('==============================================')
            self.exporting_thread = ExportingThread(
                self.ui_obj.output_directory.text(),
                self.ui_obj.export_type.currentText(),
                obj_for_export,
                self.ui_obj.progress_bar,
                self.ui_obj.progess_text)

            self.exporting_thread.prog_update.connect(self.update_progress)
            self.exporting_thread.start()
        else:
            self.ui_obj.output_directory.setStyleSheet("background: red;")

    def update_progress(self, recieved_dict):
        self.ui_obj.progress_bar.setValue(
            int((recieved_dict["Job_Number"]/recieved_dict["Total_jobs"])*100))

        if recieved_dict["Object_name"] is not None:
            self.ui_obj.progess_text.setText(f"Exporting job {recieved_dict['Job_Number']} of {recieved_dict['Total_jobs']}: {recieved_dict['Object_name']}")
        else:
            self.ui_obj.progess_text.setText("Complete!")
            # Open the directory we've just written the files to
            self.open_file_in_explorer(recieved_dict["Last_file"])

            # After X seconds close the dialog.
            countdown = 3
            for i in range(countdown, 0, -1):
                self.ui_obj.progess_text.setText(f"Closing dialog in {i-1}")
                QtWidgets.QApplication.processEvents()
                sleep(1)
            self.close_button_action()

    def browse_button_action(self, textbox_to_update):
        textbox_to_update.setStyleSheet("background: None;")
        # Get the directory location
        option = QtWidgets.QFileDialog.ShowDirsOnly
        option |= QtWidgets.QFileDialog.DontUseNativeDialog
        filename = str(QtWidgets.QFileDialog.getExistingDirectory(
            dialog, "Select Directory for output", options=option))

        if filename:
            textbox_to_update.setText(filename)

    def close_button_action(self):
        # Need to clear up the interface once a job is complete.
        self.ui_obj.output_directory.setText("")
        self.ui_obj.progress_bar.setValue(0)
        self.ui_obj.progess_text.setText("")
        for i in range(len(self.individual_obj_checkboxes)):
            self.individual_obj_checkboxes[0].deleteLater()
            del self.individual_obj_checkboxes[0]
        dialog.close()

    def open_file_in_explorer(self, file_path):
        if current_system() == 'Windows':
            Popen(r'explorer /select,"' + f'{file_path}"')
        if current_system() == 'Darwin':
            call(['open', '-R', '-aFinder', f"{file_path}/"])
