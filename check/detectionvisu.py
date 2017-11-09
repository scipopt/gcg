from Tkinter import *
from detectionplotter import *
from tkFileDialog import *
from ttk import *

class App:

    def __init__(self, master):

        self.frame = Frame(master)
        self.frame.pack()

        self.frame.winfo_toplevel().title("Plotter for detection statistics")


        for x in range(60):
            Grid.columnconfigure(self.frame, x, weight=1)

        for y in range(30):
            Grid.rowconfigure(self.frame, y, weight=1)


        self.bmaxwhite = Button(
            self.frame, text='Max white score', command=self.plotdetectionquality)
        self.bmaxwhite.grid(row=1, column=4, sticky=N+S+W+E)

        self.bmaxwhitesetpart = Button(
            self.frame, text='Max white score with setpart master', command=self.plotdetectionqualitysetpartmaster)
        self.bmaxwhitesetpart.grid(row=2, column=4, sticky=N+S+W+E)


        self.bnblocksbest = Button(
            self.frame, text='Number of blocks of whitest', command=self.plotdetectionnblocks)
        self.bnblocksbest.grid(row=1, column=3, sticky=N+S+W+E)

        self.bndecomps = Button(
            self.frame, text='Number of decompositions', command=self.plotdetectionndecomps)
        self.bndecomps.grid(row=2, column=3, sticky=N+S+W+E)

        self.bnclassesforclassifier = Button(
            self.frame, text='Number of classes for classifier', command=self.plotnclassesforclassifier)
        self.bnclassesforclassifier.grid(row=3, column=3, sticky=N+S+W+E)




    def open_file(self):
        name = askopenfilename(initialdir="/local/bastubbe/gcg-dev/check/test")
        self.plotter = Plotter(name)
        self.listboxclassifier = Listbox(self.frame)
        self.listboxclassifier.grid(row=10, column=0)

        for item in self.plotter.getclassifiernames():
            self.listboxclassifier.insert(END, item)

        self.l1text = "Filename: " + self.plotter.getfilename()
        self.l2text = "Number of instances: "+str(self.plotter.getninstances()) 
        self.l3text = "found at least one nontrivial decomp for " + str(self.plotter.getNNonTrivialDecomp()) 

        self.l4text = "found at least one nontrivial decomp with setpartitioning master for " + str(self.plotter.getNNonTrivialDecompSetpartmaster()) + "     " 


        Label(self.frame, text=self.l1text).grid(row=0, column=0, sticky=W)
        Label(self.frame, text=self.l2text).grid(row=1, column=0, sticky=W)
        Label(self.frame, text=self.l3text).grid(row=2, column=0, sticky=W)
        Label(self.frame, text=self.l4text).grid(row=3, column=0, sticky=W)



    def plotdetectionquality(self):
        self.plotter.plotdetectionquality()

    def plotdetectionqualitysetpartmaster(self):
        self.plotter.plotdetectionqualitysetpartmaster()

    def plotdetectionnblocks(self):
        self.plotter.plotnblocksofbest()

    def plotdetectionndecomps(self):
        self.plotter.plotndecomps()

    def plotnclassesforclassifier(self):
        classifierid = self.listboxclassifier.curselection()[0]
        print(classifierid)
        classifier = self.plotter.getclassifiernames()[classifierid]
        self.plotter.plotnclassesforclassifier(classifier)


root = Tk()
root.style = Style()
#('clam', 'alt', 'default', 'classic')
root.style.theme_use("clam")
menubar = Menu(root)
app = App(root)

# create a pulldown menu, and add it to the menu bar
filemenu = Menu(menubar, tearoff=0)
filemenu.add_command(label="Open", command=app.open_file)
filemenu.add_separator()
filemenu.add_command(label="Exit", command=root.quit)
menubar.add_cascade(label="File", menu=filemenu)
root.config(menu=menubar)

root.mainloop()
root.destroy() # optional; see description below
