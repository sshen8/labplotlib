# Original figure and code from https://github.com/jeffverboon/blood_expressions_plots

import matplotlib.pyplot as plt
import matplotlib.path as mpath
import matplotlib.patches as mpatches
import xml.etree.ElementTree as ET
import matplotlib.colors as colors
import os

CELL_TO_ID = {
	"GRAN": [162, 163],
	"ERY": [106, 256, 257, 254, 255, 110, 108],
	"PLT": [270, 272, 273, 274, 275, 266, 237, 235, 268, 239],
	"MDC": [157],
	"MEGA": [1],
	"CD8": [221],
	"MONO": [100],
	"pDC": [153],
	"CLP": [66],
	"HSC": [47],
	"LMPP": [53],
	"MPP": [49],
	"CMP": [68],
	"MEP": [72],
	"CD4": [78],
	"B": [88],
	"NK": [93],
	"GMPA": [122],
	"GMPB": [135],
	"GMPC": [70],
}

# reverse mapping from id to cell
ID_TO_CELL = {}
for cell, ids in CELL_TO_ID.items():
	for id in ids:
		ID_TO_CELL[id] = cell

def blood_cell_plot(data, ax=None, cmap=plt.cm.Reds, norm=None, hue=None, shapes_fname=None):
	if shapes_fname is not None:
		myshape = ET.parse(shapes_fname)
	else:
		try:
			myshape = ET.parse(os.path.join(os.path.dirname(__file__), "data", "blank_cells.eps.xml"))
		except:
			myshape = ET.parse("labplotlib/blood/data/blank_cells.eps.xml")

	if ax is None:
		fig, ax = plt.subplots()
	if norm is None:
		norm = colors.Normalize(vmin=data.min(), vmax=data.max())
	if hue is not None:
		data = data[hue]
	else:
		try:
			hue = data.name
		except:
			pass

	for path in myshape.findall("path"):
		path_data = []
		if path.attrib["type"] == "fill":
			context = dict(facecolor="none", edgecolor="black")
		elif path.attrib["type"] == "stroke":
			context = dict(facecolor="none", edgecolor="black")
		else:
			raise ValueError("Unknown type")
		id = int(path.attrib["id"])
		if id in ID_TO_CELL:
			context["facecolor"] = cmap(norm(data[ID_TO_CELL[id]]))
		for item in path:
			if item.tag == "context":
				if item.find("style").attrib["lty"] != "":
					context["linestyle"] = "dashed"
				continue
			path_data.append((
				(float(item.attrib["x"]), float(item.attrib["y"])),
				dict(
					move=mpath.Path.MOVETO,
					line=mpath.Path.LINETO,
				)[item.tag]))
		patch = mpatches.PathPatch(mpath.Path(*zip(*path_data)), **context)
		ax.add_patch(patch)
	for text in myshape.findall("text"):
		ax.text(float(text.attrib["x"]), float(text.attrib["y"]), text.attrib["string"], fontsize=10)
	ax.axis('equal') # needs to be after everything is drawn so that bounds are autoset
	ax.axis('off')

	cax = ax.inset_axes([0.07, 0.5, 0.02, 0.3]) # x0, y0, width, height
	plt.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cmap), cax=cax)
	if hue:
		cax.set_title(hue, fontsize=plt.rcParams["axes.labelsize"])
