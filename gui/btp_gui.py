#!/usr/bin/env python

import glob
import pygtk
pygtk.require('2.0')
import gtk
import ConfigParser

CONFIG_PATH = 'indurator.cfg'

class Config:
    def __init__(config_path):
        self.config = ConfigParser.SafeConfirParser()
        self.config.read(config_path)
        self.values = dict(self.config.items('main'))

    def get(var):
        if var in self.values.items():
            return self.values[var]
        else:
            return None

    def set(var, value):
        try:
            self.confg.set('main', var, value)
        except Error as e:
            raise e


class Indurator_UI(gtk.Window):
    def hello(self, widget, data=None):
        print "Hello World"

    def delete_event(self, widget, event, data=None):
        print "Delete event occured"
        return False

    def destroy(self,widget, data=None):
        gtk.main_quit()

    def toggle_snap(self, widget, spin):
        spin.set_snap_to_ticks(widget.get_active())

    def toggle_numeric(self, widget, spin):
        spin.set_numeric(widget.get_active())

    def change_digits(self, widget, spin, spin1):
        spin1.set_digits(spin.get_value_as_int())

    def get_value(self, widget, data, spin, spin2, label):
        if data == 1:
            buf = "%d" % spin.get_value_as_int()
        else:
            buf = "%0.*f" % (spin2.get_value_as_int(),
                             spin.get_value())
        label.set_text(buf)

    def reset_values(self):
        pass

    def __init__(self):
        # Initialize an empty window
        super(Indurator_UI, self).__init__()

        # Bind delete and destroy events of self
        self.connect("delete_event", self.delete_event)
        self.connect("destroy", self.destroy)

        # Set self properties
        self.set_border_width(20)
        self.set_title("Virtual Indurator")
        self.set_position(gtk.WIN_POS_CENTER_ALWAYS)

        # Add widgets to self
        main_vbox = gtk.VBox(False, 5)
        main_vbox.set_border_width(10)

        frame = gtk.Frame("Furnace Dimensions")
        main_vbox.pack_start(frame, True, True, 0)
        self.add(main_vbox)

        vbox = gtk.VBox(False, 0)
        vbox.set_border_width(5)
        frame.add(vbox)

        # Day, month, year spinners
        hbox = gtk.HBox(False, 0)
        vbox.pack_start(hbox, True, True, 5)

        vbox2 = gtk.VBox(False, 0)
        hbox.pack_start(vbox2, True, True, 5)

        label = gtk.Label("Grate Length :")
        label.set_alignment(0, 0.5)
        vbox2.pack_start(label, False, True, 0)

        adj = gtk.Adjustment(1.0, 1.0, 31.0, 1.0, 5.0, 0.0)
        spinner = gtk.SpinButton(adj, 0, 0)
        spinner.set_wrap(True)
        vbox2.pack_start(spinner, False, True, 0)

        vbox2 = gtk.VBox(False, 0)
        hbox.pack_start(vbox2, True, True, 5)

        label = gtk.Label("Grate Area :")
        label.set_alignment(0, 0.5)
        vbox2.pack_start(label, False, True, 0)

        adj = gtk.Adjustment(1.0, 1.0, 12.0, 1.0, 5.0, 0.0)
        spinner = gtk.SpinButton(adj, 0, 0)
        spinner.set_wrap(True)
        vbox2.pack_start(spinner, False, True, 0)

        hbox2 = gtk.HBox(False, 0)
        vbox.pack_start(hbox2, True, True, 5)

        vbox2 = gtk.VBox(False, 0)
        hbox2.pack_start(vbox2, True, True, 5)

        label = gtk.Label("No. of windboxes :")
        label.set_alignment(0, 0.5)
        vbox2.pack_start(label, False, True, 0)

        adj = gtk.Adjustment(1998.0, 0.0, 2100.0, 1.0, 100.0, 0.0)
        spinner = gtk.SpinButton(adj, 0, 0)
        spinner.set_wrap(False)
        spinner.set_size_request(55, -1)
        vbox2.pack_start(spinner, False, True, 0)

        frame = gtk.Frame("Pellet/Bed Properties")
        main_vbox.pack_start(frame, True, True, 0)

        vbox = gtk.VBox(False, 0)
        vbox.set_border_width(5)
        frame.add(vbox)

        hbox = gtk.HBox(False, 0)
        vbox.pack_start(hbox, False, True, 5)

        vbox2 = gtk.VBox(False, 0)
        hbox.pack_start(vbox2, True, True, 5)

        label = gtk.Label("Mean Pellet Diameter(mm) ")
        label.set_alignment(0, 0.5)
        vbox2.pack_start(label, False, True, 0)

        adj = gtk.Adjustment(0.0, -10000.0, 10000.0, 0.5, 100.0, 0.0)
        spinner1 = gtk.SpinButton(adj, 1.0, 2)
        spinner1.set_wrap(True)
        spinner1.set_size_request(100, -1)
        vbox2.pack_start(spinner1, False, True, 0)

        vbox2 = gtk.VBox(False, 0)
        hbox.pack_start(vbox2, True, True, 5)

        label = gtk.Label("Pellet Moisture Content(wt%)")
        label.set_alignment(0, 0.5)
        vbox2.pack_start(label, False, True, 0)

        adj = gtk.Adjustment(2, 1, 5, 1, 1, 0)
        spinner2 = gtk.SpinButton(adj, 0.0, 0)
        spinner2.set_wrap(True)
        vbox2.pack_start(spinner2, False, True, 0)

        hbox2 = gtk.HBox(False, 0)
        vbox.pack_start(hbox2, False, True, 5)

        vbox2 = gtk.VBox(False, 0)
        hbox2.pack_start(vbox2, True, True, 5)

        label = gtk.Label("Critical Pellet Moisture Content(wt%) ")
        label.set_alignment(0, 0.5)
        vbox2.pack_start(label, False, True, 0)

        adj = gtk.Adjustment(0.0, -10000.0, 10000.0, 0.5, 100.0, 0.0)
        spinner1 = gtk.SpinButton(adj, 1.0, 2)
        spinner1.set_wrap(True)
        spinner1.set_size_request(100, -1)
        vbox2.pack_start(spinner1, False, True, 0)

        vbox2 = gtk.VBox(False, 0)
        hbox2.pack_start(vbox2, True, True, 5)

        label = gtk.Label("Coke mean size(mm)")
        label.set_alignment(0, 0.5)
        vbox2.pack_start(label, False, True, 0)

        adj = gtk.Adjustment(2, 1, 5, 1, 1, 0)
        spinner2 = gtk.SpinButton(adj, 0.0, 0)
        spinner2.set_wrap(True)
        vbox2.pack_start(spinner2, False, True, 0)

        hbox3 = gtk.HBox(False, 0)
        vbox.pack_start(hbox3, False, True, 5)

        vbox2 = gtk.VBox(False, 0)
        hbox3.pack_start(vbox2, True, True, 5)

        label = gtk.Label("Limestone mean size(mm) ")
        label.set_alignment(0, 0.5)
        vbox2.pack_start(label, False, True, 0)

        adj = gtk.Adjustment(0.0, -10000.0, 10000.0, 0.5, 100.0, 0.0)
        spinner1 = gtk.SpinButton(adj, 1.0, 2)
        spinner1.set_wrap(True)
        spinner1.set_size_request(100, -1)
        vbox2.pack_start(spinner1, False, True, 0)

        vbox2 = gtk.VBox(False, 0)
        hbox3.pack_start(vbox2, True, True, 5)

        label = gtk.Label("Coke breeze(wt%)")
        label.set_alignment(0, 0.5)
        vbox2.pack_start(label, False, True, 0)

        adj = gtk.Adjustment(2, 1, 5, 1, 1, 0)
        spinner2 = gtk.SpinButton(adj, 0.0, 0)
        spinner2.set_wrap(True)
        vbox2.pack_start(spinner2, False, True, 0)

        frame = gtk.Frame("Process Conditions")
        main_vbox.pack_start(frame, True, True, 0)

        vbox = gtk.VBox(False, 0)
        vbox.set_border_width(5)
        frame.add(vbox)

        hbox = gtk.HBox(False, 0)
        vbox.pack_start(hbox, True, True, 5)

        vbox2 = gtk.VBox(False, 0)
        hbox.pack_start(vbox2, True, True, 5)

        label = gtk.Label("Grate Speed :")
        label.set_alignment(0, 0.5)
        vbox2.pack_start(label, False, True, 0)

        adj = gtk.Adjustment(1.0, 1.0, 31.0, 1.0, 5.0, 0.0)
        spinner = gtk.SpinButton(adj, 0, 0)
        spinner.set_wrap(True)
        vbox2.pack_start(spinner, False, True, 0)

        vbox2 = gtk.VBox(False, 0)
        hbox.pack_start(vbox2, True, True, 5)

        label = gtk.Label("Bed Height :")
        label.set_alignment(0, 0.5)
        vbox2.pack_start(label, False, True, 0)

        adj = gtk.Adjustment(1.0, 1.0, 12.0, 1.0, 5.0, 0.0)
        spinner = gtk.SpinButton(adj, 0, 0)
        spinner.set_wrap(True)
        vbox2.pack_start(spinner, False, True, 0)

        hbox2 = gtk.HBox(False, 0)
        vbox.pack_start(hbox2, True, True, 5)

        vbox2 = gtk.VBox(False, 0)
        hbox2.pack_start(vbox2, True, True, 5)

        label = gtk.Label("Hearth Layer Height :")
        label.set_alignment(0, 0.5)
        vbox2.pack_start(label, False, True, 0)

        adj = gtk.Adjustment(1998.0, 0.0, 2100.0, 1.0, 100.0, 0.0)
        spinner = gtk.SpinButton(adj, 0, 0)
        spinner.set_wrap(False)
        spinner.set_size_request(55, -1)
        vbox2.pack_start(spinner, False, True, 0)

        frame = gtk.Frame("Recuperating Hood temperature profile")
        main_vbox.pack_start(frame, True, True, 0)

        vbox = gtk.VBox(False, 0)
        vbox.set_border_width(5)
        frame.add(vbox)

        hbox = gtk.HBox(False, 0)
        vbox.pack_start(hbox, True, True, 5)

        vbox2 = gtk.VBox(False, 0)
        hbox.pack_start(vbox2, True, True, 5)

        label = gtk.Label("Hood temperature at 46m")
        label.set_alignment(0, 0.5)
        vbox2.pack_start(label, False, True, 0)

        adj = gtk.Adjustment(1.0, 1.0, 31.0, 1.0, 5.0, 0.0)
        spinner = gtk.SpinButton(adj, 0, 0)
        spinner.set_wrap(True)
        vbox2.pack_start(spinner, False, True, 0)

        vbox2 = gtk.VBox(False, 0)
        hbox.pack_start(vbox2, True, True, 5)

        label = gtk.Label("Hood temperature at 50m")
        label.set_alignment(0, 0.5)
        vbox2.pack_start(label, False, True, 0)

        adj = gtk.Adjustment(1.0, 1.0, 12.0, 1.0, 5.0, 0.0)
        spinner = gtk.SpinButton(adj, 0, 0)
        spinner.set_wrap(True)
        vbox2.pack_start(spinner, False, True, 0)

        hbox2 = gtk.HBox(False, 0)
        vbox.pack_start(hbox2, True, True, 5)

        vbox2 = gtk.VBox(False, 0)
        hbox2.pack_start(vbox2, True, True, 5)

        label = gtk.Label("Hood temperature at 100m")
        label.set_alignment(0, 0.5)
        vbox2.pack_start(label, False, True, 0)

        adj = gtk.Adjustment(1998.0, 0.0, 2100.0, 1.0, 100.0, 0.0)
        spinner = gtk.SpinButton(adj, 0, 0)
        spinner.set_wrap(False)
        spinner.set_size_request(55, -1)
        vbox2.pack_start(spinner, False, True, 0)

        vbox2 = gtk.VBox(False, 0)
        hbox2.pack_start(vbox2, True, True, 5)

        label = gtk.Label("Hood temperature at 150m")
        label.set_alignment(0, 0.5)
        vbox2.pack_start(label, False, True, 0)

        adj = gtk.Adjustment(1998.0, 0.0, 2100.0, 1.0, 100.0, 0.0)
        spinner = gtk.SpinButton(adj, 0, 0)
        spinner.set_wrap(False)
        spinner.set_size_request(55, -1)
        vbox2.pack_start(spinner, False, True, 0)


        hbox = gtk.HBox(False, 0)
        main_vbox.pack_start(hbox, False, True, 0)

        button = gtk.Button("Reset")
        button.set_size_request(80, 40)
        button.connect("clicked", self.reset_values)
        hbox.pack_start(button, True, True, 5)

        button2 = gtk.Button("Calculate")
        button2.set_size_request(80, 40)
        button2.connect("clicked", self.calculate)
        hbox.pack_start(button2, True, True, 5)
        self.show_all()

    def calculate(self, *args):
        self.hide()
        gtk.main_quit()
        Result()
        gtk.main()


class Result(gtk.Window):
    def view_graph(self, widget, event, img):
        graph = gtk.Window()
        graph.set_position(gtk.WIN_POS_CENTER_ALWAYS)
#        graph.set_size_request(400, 400)
        graph.set_title("Graph")
        vbox = gtk.VBox(False, 10)
        image = gtk.Image()
        image.set_from_file(img)
        vbox.add(image)
        graph.add(vbox)
        graph.show_all()
        return True

    def __init__(self):
        super(Result, self).__init__()
        self.connect('destroy', gtk.main_quit)

        self.set_border_width(20)
        self.set_title("Result")
        self.set_position(gtk.WIN_POS_CENTER_ALWAYS)

        vbox = gtk.VBox(False, 10)

        graphs = []
        for each in glob.glob('graphs/*'):
            event_box = gtk.EventBox()
            event_box.connect("button_release_event", self.view_graph, each)
            vbox.add(event_box)
            image = gtk.Image()
            pixbuf = gtk.gdk.pixbuf_new_from_file(each)
            scaled_buf = pixbuf.scale_simple(200, 200, gtk.gdk.INTERP_BILINEAR)
            image.set_from_pixbuf(scaled_buf)
            event_box.add(image)

        self.add(vbox)
        self.show_all()

class Welcome(gtk.Window):

    def start_computation(self, *args):
        self.hide()
        gtk.main_quit()
        Indurator_UI()
        gtk.main()

    def __init__(self):
        super(Welcome, self).__init__()

        # Bind delete and destroy events of self
        self.connect("destroy", gtk.main_quit)

        # Set self properties
        self.set_border_width(10)
        self.set_title("Virtual Indurator")
        self.set_position(gtk.WIN_POS_CENTER_ALWAYS)

        vbox = gtk.VBox(False, 0)
        hbox = gtk.HBox(True, 3)

        valign = gtk.Alignment(0, 1, 0, 0)
        vbox.pack_start(valign)

        image = gtk.Image()
        image.set_from_file('indurator.png')
        image.show()
        vbox.pack_start(image, True, True, 0)

        button = gtk.Button("Start")
        button.set_size_request(80,40)
        button.connect("clicked", self.start_computation)
        button.set_tooltip_text("Click to start the Simulation")

        hbox.add(button)
        halign = gtk.Alignment(1, 0, 0, 0)
        halign.add(hbox)
        vbox.pack_start(halign, False, False, 3)
        self.add(vbox)
        self.show_all()

if __name__ == "__main__":

    Welcome()
    gtk.main()
