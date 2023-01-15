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
surface_export_types = ["Surface", "Cartoon", "Sticks", "Mesh", "Lines",
                        "Ribbon", "Wire", "Licorice", "Dots", "Spheres"]

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
# Exporting functions Thread:
# =============================================================================
class Exporting_Thread(QtCore.QThread):
    progress_bar_updating = QtCore.pyqtSignal(dict)

    def __init__(self, directory_output_location, export_type, output_obj_list, progress_bar_obj, progress_label):
        QtCore.QThread.__init__(self, None)
        self.directory_output_location = directory_output_location
        self.export_type = export_type
        self.output_obj_list = output_obj_list
        self.progress_bar_obj = progress_bar_obj
        self.progress_label = progress_label

    def run(self):
        # Hide everything in the scene so it can be individually displayed.
        cmd.hide('everything')
        print('==============================================')
        print('Creating surfaces and exporting...:')
        print('==============================================')

        # Show each object as determined by the user.
        for object_no in range(len(self.output_obj_list)):
            cmd.disable('all')
            object_name = self.output_obj_list[object_no].text()
            print(f"    -{object_name}")
            # Send over the data to the progress bar to reflect the current state
            self.progress_bar_updating.emit({"Job_Number": object_no,
                                             "Total_jobs": len(self.output_obj_list),
                                             "Object_name": object_name})
            cmd.show_as(self.export_type.lower(), object_name)
            cmd.enable(object_name)
            name_for_output = path.join(self.directory_output_location, f"{object_name}.wrl")
            cmd.save(name_for_output, ['wrl'])
            self.output_obj_list[object_no].setStyleSheet("background-color: green")
        print(' ')
        self.progress_bar_updating.emit({"Job_Number": len(self.output_obj_list),
                                         "Total_jobs": len(self.output_obj_list),
                                         "Object_name": None,
                                         "Last_file": name_for_output})

        print('==============================================')
        cmd.disable('all')
        print('Complete')
        print('==============================================')

# =============================================================================
# GUI elements:
# =============================================================================
class MainWindow(QtWidgets.QMainWindow):
    def __init__(self, *args, **kwargs):
        super(MainWindow, self).__init__(*args, **kwargs)
        uifile = path.join(path.dirname(__file__), 'Main_menu.ui')
        self.user_interface = loadUi(uifile, dialog)
        self.scroll_widget_area = QtWidgets.QWidget()
        self.v_layout = QtWidgets.QVBoxLayout()
        self.scroll_widget_area.setLayout(self.v_layout)
        self.user_interface.scrollArea.setWidget(self.scroll_widget_area)

        self.individual_obj_checkboxes = []
        self.connect_interface_buttons()
        self.populate_individual_objects()

    def connect_interface_buttons(self):
        self.user_interface.refresh_button.clicked.connect(self.populate_individual_objects)
        self.user_interface.select_all.clicked.connect(lambda: self.all_selection(select_all=True))
        self.user_interface.select_none.clicked.connect(lambda: self.all_selection(select_all=False))
        self.user_interface.browse_button.clicked.connect(lambda: self.browse_button_function(self.user_interface.output_directory))
        self.user_interface.Submit_button.clicked.connect(self.submit_button_function)
        self.user_interface.Cancel_button.clicked.connect(self.close_button_function)
        self.populate_export_types(self.user_interface.export_type)

    # ==========================================================
    # Return the dialog for handling in PyMOL
    # ==========================================================
    def return_dialog(self):
        return self.user_interface

    def populate_export_types(self, dropdown_box_to_update):
        for each_surface in surface_export_types:
            dropdown_box_to_update.addItem(each_surface)

    def all_selection(self, select_opt):
        for each_checkbox in self.individual_obj_checkboxes:
            each_checkbox.setChecked(select_opt)

    def populate_individual_objects(self):
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

    def submit_button_function(self):
        if self.user_interface.output_directory.text() != "":
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
                print(f"{len(obj_for_export)} objects selected, exporting all to {self.user_interface.output_directory.text()}")
                print('==============================================')
            self.exporting_thread = Exporting_Thread(self.user_interface.output_directory.text(),
                                                     self.user_interface.export_type.currentText(),
                                                     obj_for_export,
                                                     self.user_interface.progress_bar,
                                                     self.user_interface.progess_text)
            self.exporting_thread.progress_bar_updating.connect(self.update_progress)
            self.exporting_thread.start()
        else:
            self.user_interface.output_directory.setStyleSheet("background: red;")

    def update_progress(self, recieved_dict):
        self.user_interface.progress_bar.setValue(int((recieved_dict["Job_Number"]/recieved_dict["Total_jobs"])*100))
        if recieved_dict["Object_name"] is not None:
            self.user_interface.progess_text.setText(f"Exporting job {recieved_dict['Job_Number']} of {recieved_dict['Total_jobs']}: {recieved_dict['Object_name']}")
        else:
            self.user_interface.progess_text.setText("Complete!")
            # Open the directory we've just written the files to
            self.open_file_in_explorer(recieved_dict["Last_file"])

            # After X seconds close the dialog.
            countdown = 3
            for i in range(countdown, 0, -1):
                self.user_interface.progess_text.setText(f"Closing dialog in {i-1}")
                QtWidgets.QApplication.processEvents()
                sleep(1)
            self.close_button_function()

    def browse_button_function(self, textbox_for_updating):
        textbox_for_updating.setStyleSheet("background: None;")
        # Get the directory location
        option = QtWidgets.QFileDialog.ShowDirsOnly
        option |= QtWidgets.QFileDialog.DontUseNativeDialog
        filename = str(QtWidgets.QFileDialog.getExistingDirectory(dialog, "Select Directory for output", options=option))
        if filename:
            textbox_for_updating.setText(filename)

    def close_button_function(self):
        # Need to clear up the interface once a job is complete.
        self.user_interface.output_directory.setText("")
        self.user_interface.progress_bar.setValue(0)
        self.user_interface.progess_text.setText("")
        for i in range(len(self.individual_obj_checkboxes)):
            self.individual_obj_checkboxes[0].deleteLater()
            del self.individual_obj_checkboxes[0]
        dialog.close()

    def open_file_in_explorer(self, file_path):
        if current_system() == 'Windows':
            Popen(r'explorer /select,"' + f'{file_path}"')
        if current_system() == 'Darwin':
            call(['open', '-R', '-aFinder', f"{file_path}/"])
