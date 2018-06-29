//==============================================================================
//	
//	Copyright (c) 2002-
//	Authors:
//	* Alistair John Strachan <alistair@devzero.co.uk> (University of Edinburgh)
//	* Mike Arthur <mike@mikearthur.co.uk> (University of Edinburgh)
//	* Zak Cohen <zakcohen@gmail.com> (University of Edinburgh)
//	* Mark Kattenbelt <mark.kattenbelt@comlab.ox.ac.uk> (University of Oxford, formerly University of Birmingham)
//	* Clemens Wiltsche <clemens.wiltsche@cs.ox.ac.uk> (University of Oxford)
//	
//------------------------------------------------------------------------------
//	
//	This file is part of PRISM.
//	
//	PRISM is free software; you can redistribute it and/or modify
//	it under the terms of the GNU General Public License as published by
//	the Free Software Foundation; either version 2 of the License, or
//	(at your option) any later version.
//	
//	PRISM is distributed in the hope that it will be useful,
//	but WITHOUT ANY WARRANTY; without even the implied warranty of
//	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//	GNU General Public License for more details.
//	
//	You should have received a copy of the GNU General Public License
//	along with PRISM; if not, write to the Free Software Foundation,
//	Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
//	
//==============================================================================

package userinterface.graph;

import java.awt.*;
import java.awt.print.*;
import java.io.*;
import java.util.*;
import java.awt.image.*;
import java.awt.BasicStroke;

import javax.swing.JLabel;
import javax.swing.JSlider;
import javax.swing.JTextField;
import javax.swing.JPanel;
import javax.swing.JButton;
import javax.swing.JRadioButton;
import javax.swing.JCheckBox;
import javax.swing.Box;
import javax.swing.BoxLayout;
import javax.swing.ButtonGroup;
import javax.swing.AbstractButton;
import javax.swing.event.ChangeListener;
import javax.swing.event.ChangeEvent;

import java.awt.event.ItemListener;
import java.awt.event.ItemEvent;
import java.awt.event.ActionListener;
import java.awt.event.ActionEvent;

import javax.swing.SwingUtilities;
import javax.swing.JScrollPane;

import org.w3c.dom.DOMException;
import org.w3c.dom.Document;
import org.w3c.dom.Element;
import org.w3c.dom.NodeList;
import org.xml.sax.*;

import javax.xml.parsers.*;
import javax.xml.transform.*;
import javax.xml.transform.dom.*;
import javax.xml.transform.stream.*;

import org.jfree.ui.*;
import org.jfree.chart.*;
import org.jfree.chart.plot.*;
import org.jfree.chart.title.*;
import org.jfree.chart.axis.*;
import org.jfree.data.xy.*;
import org.jfree.chart.renderer.xy.*;
import org.jfree.chart.encoders.*;
import org.jfree.chart.event.ChartProgressEvent;
import org.jfree.chart.annotations.XYPointerAnnotation;
import org.jfree.chart.annotations.XYPolygonAnnotation;

import net.sf.epsgraphics.*;
import prism.*;
import userinterface.GUIPlugin;
import settings.*;

/**
 * This class represents multiple series on a single unit graph; all series are
 * rendered together by MultiGraphView.
 */
public class Graph extends ChartPanel implements SettingOwner, EntityResolver, Observer, Printable
{
	/** Interval for graph updates (in ms). */
	private static final long updateInterval = 250;

	/** Actual JFreeChart representation of datasets. */
	private JFreeChart chart;

	/** XYPlot of this JFreeChart */
	private XYPlot plot;

	/**
	 * List of currently available series data to plot. (Make sure to
	 * synchronize)
	 */
	private XYSeriesCollection seriesCollection;

	/**
	 * Maps SeriesKeys to a XYSeries. (Make sure to synchronize on
	 * seriesCollection)
	 */
	private HashMap<SeriesKey, XYSeries> keyToSeries;

	/**
	 * Maps SeriesKeys to a PointList. (Make sure to synchronize on
	 * seriesCollection)
	 */
	private HashMap<SeriesKey, PointList> keyToPointList;

	/**
	 * Maps SeriesKeys to a Graph Series. (Make sure to synchronize on
	 * seriesCollection)
	 */
	private HashMap<SeriesKey, SeriesSettings> keyToGraphSeries;

	/**
	 * Allows us to batch graph points (JFreeChart is not realtime). (Make sure
	 * to synchronize on seriesCollection)
	 */
	private HashMap<SeriesKey, LinkedList<XYDataItem>> graphCache;
	private HashMap<SeriesKey, LinkedList<double[]>> polyCache;

	/** Display for settings. Required to implement SettingsOwner */
	private SettingDisplay display;

	/** Settings of this graph. */
	private MultipleLineStringSetting graphTitle;
	private FontColorSetting titleFont;
	private BooleanSetting legendVisible;
	private ChoiceSetting legendPosition;
	private FontColorSetting legendFont;

	/** Settings of the axis. */
	private AxisSettings xAxisSettings;

	private AxisSettings yAxisSettings;

	/** Display settings */
	private DisplaySettings displaySettings;

	/** GraphSeriesList */
	private SeriesSettingsList seriesList;

	/** legend position */
	public static final int LEFT = 0;
	public static final int RIGHT = 1;
	public static final int BOTTOM = 2;
	public static final int TOP = 3;

	// sliders
	private JScrollPane sliderPane = null;
	private JPanel sliders = new JPanel();

	
	public JScrollPane getSliderPane()
	{
		return sliderPane;
	}

	/**
	 * Initialises the GraphModel's series and canvas list. Also starts off the
	 * graph update timer (one per chart).
	 */
	public Graph()
	{
		this("");
	}

	/**
	* Generate a graph displaying a Pareto set provided as a PointList
	**/
	public void addPointListSeries(PointList pointList)
	{
		if(keyToPointList == null) // none added yet
			addFirstPointListSeries(pointList);
		else
			addAdditionalPointListSeries(pointList);
	}
	private void addFirstPointListSeries(PointList pointList)
	{
		try {
			legendVisible.setValue(new Boolean(false)); // hide legend
			this.notifySettingChanged(legendVisible);
		} catch (SettingException e) {
			// should not happen, otherwise, ignore
		}
		XYItemRenderer r = plot.getRenderer();
		// if possible, try to match the old grapher
		if (r instanceof XYLineAndShapeRenderer) {
			XYLineAndShapeRenderer renderer = (XYLineAndShapeRenderer) r;
				renderer.setBaseShapesVisible(false);
				renderer.setBaseShapesFilled(false);
				renderer.setLinesVisible(false); // no lines connecting shown
				renderer.setShapesVisible(false); // no points shown
				renderer.setDrawSeriesLineAsPath(false);
				renderer.setAutoPopulateSeriesPaint(true);
				renderer.setAutoPopulateSeriesShape(true);
		}

		// remove all series
		synchronized (seriesCollection) {
			for(SeriesKey sk : graphCache.keySet()) {
				removeSeries(sk);
			}
		}
		
		keyToPointList = new HashMap<SeriesKey, PointList>();

		
		addAdditionalPointListSeries(pointList);

		// slice, initially picks dimensions 0 and 1, so fill in the midpoints in the remaining dimensions
		final Slice slice = new Slice();
		for (int i = 2; i < pointList.getFullDim(); i++)
			slice.put(i, pointList.getMidpoint(i));

		// make initial drawing
		redraw(slice);

		// make container for sliders
		GridBagLayout gridbag = new GridBagLayout();
		GridBagConstraints c = new GridBagConstraints();

		// set layout
		sliders.setLayout(gridbag);
		int width = 7; // layout width

		// add axis switch button
		c.fill = GridBagConstraints.NONE;
		c.gridwidth = 2;
		c.gridx = width - 2;
		c.gridy = 0;
		JButton switch_axes = new JButton("<->");
		switch_axes.setToolTipText("Switch axes");
		gridbag.setConstraints(switch_axes, c);
		sliders.add(switch_axes);

		// label for selection
		c.fill = GridBagConstraints.HORIZONTAL;
		c.weightx = 5;
		c.gridwidth = 5;
		c.gridx = 1;
		c.gridy = 1;
		JLabel slabel = new JLabel("Slice at");
		slabel.setToolTipText("Value at which Pareto set is sliced.");
		gridbag.setConstraints(slabel, c);
		sliders.add(slabel); // label the expressions
		c.fill = GridBagConstraints.NONE;

		// label for x axis
		c.weightx = 0.1;
		c.gridwidth = 1;
		c.gridx = width - 2;
		c.gridy = 1;
		JLabel xlabel = new JLabel("X");
		xlabel.setToolTipText("Select expression on abscissa.");
		gridbag.setConstraints(xlabel, c);
		sliders.add(xlabel); // label the button
		// label for y axis
		c.weightx = 0.1;
		c.gridwidth = 1;
		c.gridx = width - 1;
		c.gridy = 1;
		JLabel ylabel = new JLabel("Y");
		ylabel.setToolTipText("Select expression on ordinate.");
		gridbag.setConstraints(ylabel, c);
		sliders.add(ylabel); // label the button

		// label for projection
		c.weightx = 0.1;
		c.gridwidth = 1;
		c.gridx = width;
		c.gridy = 1;
		JLabel plabel = new JLabel("Project");
		plabel.setToolTipText("Select expressions to project.");
		gridbag.setConstraints(plabel, c);
		sliders.add(plabel); // label the box

		c.fill = GridBagConstraints.HORIZONTAL;


		final java.util.List<JTextField> sliceSelect = new ArrayList<JTextField>();
		final java.util.List<JCheckBox> projectSelect = new ArrayList<JCheckBox>();
		ButtonGroup xAxisSelect = new ButtonGroup(); // radio buttons for selecting x axis
		ButtonGroup yAxisSelect = new ButtonGroup(); // radio buttons for selecting y axis

		final java.util.List<JRadioButton> xbuttons = new ArrayList<JRadioButton>();
		final java.util.List<JRadioButton> ybuttons = new ArrayList<JRadioButton>();

		switch_axes.addActionListener(new ActionListener()
		{
			@Override
			public void actionPerformed(ActionEvent e)
			{
				// get selected xbutton
				int ix = 0;
				for (ix = 0; ix < xbuttons.size(); ix++)
					if (xbuttons.get(ix).isSelected())
						break;

				// get selected ybutton
				int iy = 0;
				for (iy = 0; ix < ybuttons.size(); iy++)
					if (ybuttons.get(iy).isSelected())
						break;

				// change around now
				xbuttons.get(iy).setEnabled(true);
				xbuttons.get(iy).setSelected(true);
				xbuttons.get(ix).setEnabled(false);
				ybuttons.get(ix).setEnabled(true);
				ybuttons.get(ix).setSelected(true);
				ybuttons.get(iy).setEnabled(false);

				// need to wait for sliders to be enabled/disabled
				SwingUtilities.invokeLater(new Runnable()
				{
					public void run()
					{
						// recalculate the slice
						reSlice(slice, sliceSelect, projectSelect, xbuttons, ybuttons);
						redraw(slice);
					}
				});
			}
		});

		for (int j = 0; j < pointList.getFullDim(); j++) { // add options to choose slice
			final int slider_i = j;
			final JTextField manual = new JTextField(); // textbox for manual editing
			sliceSelect.add(manual);
			final JCheckBox projectbox = new JCheckBox(); // check box to project out dimensions
			projectSelect.add(projectbox);

			// radio buttons for selecting sliders
			final JRadioButton xbutton = new JRadioButton();
			final JRadioButton ybutton = new JRadioButton();
			xbuttons.add(xbutton);
			ybuttons.add(ybutton);
			xbutton.addChangeListener(new ChangeListener()
			{
				@Override
				public void stateChanged(ChangeEvent e)
				{
					if (xbutton.isSelected()) {
						ybutton.setEnabled(false); // if selected as x, cannot select as y
						manual.setEnabled(false);
						projectbox.setEnabled(false);
						projectbox.setSelected(false);
					} else {
						ybutton.setEnabled(true); // if no longer selected as x, can select as y
						if (ybutton.isSelected()) {
							manual.setEnabled(false);
							projectbox.setEnabled(false);
							projectbox.setSelected(false);
						} else {
							manual.setEnabled(true);
							projectbox.setEnabled(true);
						}
					}
					// need to wait for sliders to be enabled/disabled
					SwingUtilities.invokeLater(new Runnable()
					{
						public void run()
						{
							// recalculate the slice
							reSlice(slice, sliceSelect, projectSelect, xbuttons, ybuttons);
							redraw(slice);
						}
					});
				}
			});
			ybutton.addChangeListener(new ChangeListener()
			{
				@Override
				public void stateChanged(ChangeEvent e)
				{
					if (ybutton.isSelected()) {
						xbutton.setEnabled(false); // if selected as y, cannot select as x
						manual.setEnabled(false);
						projectbox.setEnabled(false);
						projectbox.setSelected(false);
					} else {
						xbutton.setEnabled(true); // if no longer selected as y, can select as x
						if (xbutton.isSelected()) {
							manual.setEnabled(false);
							projectbox.setEnabled(false);
							projectbox.setSelected(false);
						} else {
							manual.setEnabled(true);
							projectbox.setEnabled(true);
						}
					}
					// need to wait for sliders to be enabled/disabled
					SwingUtilities.invokeLater(new Runnable()
					{
						public void run()
						{
							// recalculate the slice
							reSlice(slice, sliceSelect, projectSelect, xbuttons, ybuttons);
							redraw(slice);
						}
					});
				}
			});

			// initial defaults
			if (j == 0)
				xbutton.setSelected(true);
			if (j == 1)
				ybutton.setSelected(true);

			// group the buttons
			xAxisSelect.add(xbutton); // add to buttons selecting x axis
			yAxisSelect.add(ybutton); // add to buttons selecting y axis

			// put slider in layout, with label and text field
			c.weightx = 0.5;
			c.gridwidth = 1;
			c.gridx = 0;
			c.gridy = j + 2;
			JLabel label = new JLabel(pointList.getLabel(slider_i), JLabel.CENTER);
			gridbag.setConstraints(label, c);
			sliders.add(label); // label the slider

			c.weightx = 5;
			c.gridwidth = width - 3;
			c.gridx = 1;
			c.gridy = j + 2;

			gridbag.setConstraints(manual, c);
			sliders.add(manual); // add manual edit text field
			manual.setToolTipText("Manually enter slice position for " + pointList.getLabel(slider_i) + ".");
			manual.setText(String.format("%.2f", pointList.getMidpoint(slider_i)));
			manual.addActionListener(new ActionListener()
			{
				public void actionPerformed(ActionEvent e)
				{
					try {
						// get average midpoint of all point lists
						// and the min min and max max of all point lists
						double val = 0.0, minmin = 0.0, maxmax = 0.0;
						
						synchronized (seriesCollection) {
							for (SeriesKey sk : graphCache.keySet()) {
								PointList pointList = keyToPointList.get(sk);
								val += pointList.getMidpoint(slider_i);
								minmin = minmin <= pointList.getCutoffMin(slider_i) ? minmin : pointList.getCutoffMin(slider_i);
								maxmax = maxmax >= pointList.getCutoffMax(slider_i) ? maxmax : pointList.getCutoffMax(slider_i);
							}
							val /= graphCache.size(); // take average
						}
						try {
							val = Double.parseDouble(manual.getText());
						} catch (NumberFormatException ey) {
							manual.setText(String.format("%.2f", val));
						}
						if (val <= minmin) {
							val = minmin;
							manual.setText(String.format("%.2f", val));
						} else if (val >= maxmax) {
							val = minmin;
							manual.setText(String.format("%.2f", val));
						}
						reSlice(slice, sliceSelect, projectSelect, xbuttons, ybuttons);
						redraw(slice);
					} catch (NumberFormatException ex) {
						// ignore
					}
				}
			});

			// add radio buttons to layout
			c.fill = GridBagConstraints.NONE;
			c.weightx = 0.1;
			c.gridwidth = 1;
			c.gridx = width - 2;
			c.gridy = j + 2;
			xbutton.setToolTipText("Set " + pointList.getLabel(slider_i) + " to be on the abscissa.");
			gridbag.setConstraints(xbutton, c);
			sliders.add(xbutton);
			ybutton.setToolTipText("Set " + pointList.getLabel(slider_i) + " to be on the ordinate.");
			c.weightx = 0.1;
			c.gridwidth = 1;
			c.gridx = width - 1;
			c.gridy = j + 2;
			gridbag.setConstraints(ybutton, c);
			sliders.add(ybutton);

			// add projection box to layout
			c.weightx = 1;
			c.gridwidth = 1;
			c.gridx = width;
			c.gridy = j + 2;
			gridbag.setConstraints(projectbox, c);
			projectbox.setToolTipText("Project on " + pointList.getLabel(slider_i) + ".");
			sliders.add(projectbox);
			projectbox.addItemListener(new ItemListener()
			{
				public void itemStateChanged(ItemEvent e)
				{
					// first adjust enabled status
					if (projectbox.isSelected()) {
						if (xbutton.isSelected() || ybutton.isSelected()) { // sanity check throws SettingsException
							projectbox.setSelected(false); // don't project out selected dimensions
							return;
						}
						xbutton.setEnabled(false);
						ybutton.setEnabled(false);
						manual.setEnabled(false);
					} else {
						if (!xbutton.isSelected() && !ybutton.isSelected()) { // sanity check
							xbutton.setEnabled(true);
							ybutton.setEnabled(true);
							manual.setEnabled(true);
						}
					}
					// now reslice and reproject
					reSlice(slice, sliceSelect, projectSelect, xbuttons, ybuttons);
					redraw(slice);
				}
			});
			c.fill = GridBagConstraints.HORIZONTAL;
		}
		
		this.sliderPane.revalidate();
		this.sliderPane.repaint();
	}

	// recomputes the slice after changing the axes
	private void reSlice(Slice slice, java.util.List<JTextField> sliceSelect, java.util.List<JCheckBox> projectSelect, java.util.List<JRadioButton> xbuttons,
			java.util.List<JRadioButton> ybuttons)
	{
		slice.clear();

		int ix = 0;
		int iy = 0;

		// fill the slice with the values of the sliders
		for (int i = 0; i < sliceSelect.size(); i++) {
			JTextField manual = sliceSelect.get(i);
			if (manual.isEnabled()) {
				// get average midpoint of all point lists
				double val = 0.0;				
				synchronized (seriesCollection) {
					for (SeriesKey sk : graphCache.keySet()) {
						PointList pointList = keyToPointList.get(sk);
						val += pointList.getMidpoint(i);
					}
					val /= graphCache.size(); // take average
				}
				try {
					val = Double.parseDouble(manual.getText());
				} catch (NumberFormatException e) {
					manual.setText(String.format("%.2f", val));
				}
				slice.put(i, val);
			} else {
				if (xbuttons.get(i).isSelected())
					ix = i;
				if (ybuttons.get(i).isSelected())
					iy = i;
				if (projectSelect.get(i).isSelected())
					slice.project(i);
			}
		}

		if (ix > iy)
			slice.swapdimensions = true;
	}

	private void redraw(Slice slice)
	{
		synchronized (seriesCollection) {
			for (SeriesKey sk : graphCache.keySet()) {
				PointList pointList = keyToPointList.get(sk);

				// remove all data from series for now:
				removeAllPoints(sk);

				// update pointList with the given slice
				try {
					pointList.updateParetoSet(slice);
				} catch (PrismException exp) {
					exp.printStackTrace();
				}

				//Get points in pointlist --- they form the bounding rectangle in the display
				java.util.List<prism.Point> l = pointList.getPoints();
				for (prism.Point p : l)
					addPointToSeries(sk, new XYDataItem(p.getCoord(0), p.getCoord(1)));

				// redraw polyhedra as annotations
				plot.clearAnnotations();
				java.util.List<java.util.List<prism.Point>> polys = pointList.getPolyPoints();
				for (java.util.List<prism.Point> poly : polys) {
					double[] polyReal = new double[poly.size() * 2];
					for (int i = 0; i < poly.size(); i++) {
						prism.Point p = poly.get(i);
						polyReal[2 * i] = p.getCoord(0);
						polyReal[2 * i + 1] = p.getCoord(1);
					}
					addPolyToSeries(sk, polyReal);
				}

				// set axis labels
				getXAxisSettings().setHeading(pointList.getXLabel());
				getYAxisSettings().setHeading(pointList.getYLabel());

			}
		}
	}

	
	/**
	 * Initialises the GraphModel's series and canvas list. Also starts off the
	 * graph update timer (one per chart).
	 * 
	 * @param title
	 *            Title of the graph.
	 */
	public Graph(String title)
	{
		super(ChartFactory.createXYLineChart(title, "X", "Y", new XYSeriesCollection(), PlotOrientation.VERTICAL, true, true, false));

		// add sliders to graph panel - invisible for now
		sliderPane = new JScrollPane(sliders);
		
		graphCache = new HashMap<SeriesKey, LinkedList<XYDataItem>>();
		polyCache = new HashMap<SeriesKey, LinkedList<double[]>>();
		keyToSeries = new HashMap<SeriesKey, XYSeries>();
		keyToGraphSeries = new HashMap<SeriesKey, SeriesSettings>();
		graphTitle = new MultipleLineStringSetting("title", title, "The main title heading for the chart.", this, false);
		titleFont = new FontColorSetting("title font", new FontColorPair(new Font("SansSerif", Font.PLAIN, 14), Color.black), "The font for the chart's title",
				this, false);
		legendVisible = new BooleanSetting("legend visible?", new Boolean(true),
				"Should the legend, which displays all of the series headings, be displayed?", this, false);

		String[] choices = { "Left", "Right", "Bottom", "Top" };
		legendPosition = new ChoiceSetting("legend position", choices, choices[RIGHT], "The position of the legend", this, false);
		legendFont = new FontColorSetting("legend font", new FontColorPair(new Font("SansSerif", Font.PLAIN, 11), Color.black), "The font for the legend",
				this, false);

		// Some easy references
		chart = super.getChart();
		plot = chart.getXYPlot();
		plot.setBackgroundPaint((Paint) Color.white);
		seriesCollection = (XYSeriesCollection) plot.getDataset();

		xAxisSettings = new AxisSettings("X", true, this);
		yAxisSettings = new AxisSettings("Y", false, this);

		xAxisSettings.addObserver(this);
		yAxisSettings.addObserver(this);

		displaySettings = new DisplaySettings(this);
		displaySettings.addObserver(this);

		seriesList = new SeriesSettingsList(this);

		// create a regular XY line chart
		XYItemRenderer r = plot.getRenderer();
		// if possible, try to match the old grapher
		if (r instanceof XYLineAndShapeRenderer) {
			XYLineAndShapeRenderer renderer = (XYLineAndShapeRenderer) r;
				renderer.setBaseShapesVisible(true);
				renderer.setBaseShapesFilled(true);
				renderer.setAutoPopulateSeriesPaint(true);
				renderer.setAutoPopulateSeriesShape(true);
		}

		plot.setDrawingSupplier(new DefaultDrawingSupplier(SeriesSettings.DEFAULT_PAINTS, DefaultDrawingSupplier.DEFAULT_OUTLINE_PAINT_SEQUENCE,
				DefaultDrawingSupplier.DEFAULT_STROKE_SEQUENCE, DefaultDrawingSupplier.DEFAULT_OUTLINE_STROKE_SEQUENCE, SeriesSettings.DEFAULT_SHAPES));

		super.setPopupMenu(null);

		/* Make sure the graph resembles its default settings. */
		updateGraph();

		// schedule a periodic timer for graph updates
		new java.util.Timer().scheduleAtFixedRate(new GraphUpdateTask(), 0, // start now
				updateInterval);
	}

	public int compareTo(Object o)
	{
		if (o instanceof SettingOwner) {
			SettingOwner po = (SettingOwner) o;
			if (getSettingOwnerID() < po.getSettingOwnerID())
				return -1;
			else if (getSettingOwnerID() > po.getSettingOwnerID())
				return 1;
			else
				return 0;
		} else
			return 0;
	}

	public int getNumSettings()
	{
		return 5;
	}

	public Setting getSetting(int index)
	{
		switch (index) {
		case 0:
			return graphTitle;
		case 1:
			return titleFont;
		case 2:
			return legendVisible;
		case 3:
			return legendPosition;
		case 4:
			return legendFont;
		default:
			return null;
		}
	}

	public String getSettingOwnerClassName()
	{
		return "Model";
	}

	public int getSettingOwnerID()
	{
		return prism.PropertyConstants.MODEL;
	}

	public String getSettingOwnerName()
	{
		return graphTitle.getStringValue();
	}

	public void doEnables()
	{
		legendPosition.setEnabled(legendVisible.getBooleanValue());
		legendFont.setEnabled(legendVisible.getBooleanValue());
	}

	public void update(Observable o, Object arg)
	{
		if (o == xAxisSettings) {
			/* X axis changed */
			super.repaint();
		} else if (o == yAxisSettings) {
			/* Y axis changed */
			super.repaint();
		} else if (o == displaySettings) {
			/* Display settings changed */
			super.repaint();
		} else {
			for (Map.Entry<SeriesKey, SeriesSettings> entry : keyToGraphSeries.entrySet()) {
				/* Graph series settings changed */
				if (entry.getValue().equals(o))
					repaint();
			}
		}
	}

	public void notifySettingChanged(Setting setting)
	{
		updateGraph();
	}

	private void updateGraph()
	{
		/* Update title if necessary. */
		if (!this.chart.getTitle().equals(graphTitle)) {
			this.chart.setTitle(graphTitle.getStringValue());
		}

		/* Update title font if necessary. */
		if (!titleFont.getFontColorValue().f.equals(this.chart.getTitle().getFont())) {
			this.chart.getTitle().setFont(titleFont.getFontColorValue().f);
		}

		/* Update title colour if necessary. */
		if (!titleFont.getFontColorValue().c.equals(this.chart.getTitle().getPaint())) {
			this.chart.getTitle().setPaint(titleFont.getFontColorValue().c);
		}

		if (legendVisible.getBooleanValue() != (this.chart.getLegend() != null)) {
			if (!legendVisible.getBooleanValue()) {
				this.chart.removeLegend();
			} else {
				LegendTitle legend = new LegendTitle(plot.getRenderer());
				legend.setBackgroundPaint(Color.white);
				legend.setBorder(1, 1, 1, 1);

				this.chart.addLegend(legend);
			}
		}

		if (this.chart.getLegend() != null) {
			LegendTitle legend = this.chart.getLegend();

			/* Put legend on the left if appropriate. */
			if ((legendPosition.getCurrentIndex() == LEFT) && !legend.getPosition().equals(RectangleEdge.LEFT)) {
				legend.setPosition(RectangleEdge.LEFT);
			}
			/* Put legend on the right if appropriate. */
			if ((legendPosition.getCurrentIndex() == RIGHT) && !legend.getPosition().equals(RectangleEdge.RIGHT)) {
				legend.setPosition(RectangleEdge.RIGHT);
			}
			/* Put legend on the top if appropriate. */
			if ((legendPosition.getCurrentIndex() == TOP) && !legend.getPosition().equals(RectangleEdge.TOP)) {
				legend.setPosition(RectangleEdge.TOP);
			}
			/* Put legend on the bottom if appropriate. */
			if ((legendPosition.getCurrentIndex() == BOTTOM) && !legend.getPosition().equals(RectangleEdge.BOTTOM)) {
				legend.setPosition(RectangleEdge.BOTTOM);
			}

			/* Set legend font. */
			if (!legend.getItemFont().equals(legendFont.getFontColorValue().f)) {
				legend.setItemFont(legendFont.getFontColorValue().f);
			}
			/* Set legend font colour. */
			if (!legend.getItemPaint().equals(legendFont.getFontColorValue().c)) {
				legend.setItemPaint(legendFont.getFontColorValue().c);
			}
		}

		super.repaint();
		doEnables();
	}

	public void setDisplay(SettingDisplay display)
	{
		this.display = display;
	}

	public SettingDisplay getDisplay()
	{
		return display;
	}

	/** 
	 * Returns an object that you have to synchronise in one case:
	 * - You depend on series not changing.
	 */
	public Object getSeriesLock()
	{
		return seriesCollection;
	}

	public java.util.Vector<SeriesKey> getAllSeriesKeys()
	{
		synchronized (seriesCollection) {
			java.util.Vector<SeriesKey> result = new java.util.Vector<SeriesKey>();

			for (Map.Entry<SeriesKey, XYSeries> entries : keyToSeries.entrySet()) {
				result.add(entries.getKey());
			}

			return result;
		}
	}

	public SeriesSettingsList getGraphSeriesList()
	{
		return seriesList;
	}

	/**
	 * Should always be synchronised on seriesCollection when called.
	 */
	public SeriesSettings getGraphSeries(SeriesKey key)
	{
		synchronized (seriesCollection) {
			if (keyToGraphSeries.containsKey(key)) {
				return keyToGraphSeries.get(key);
			}

			return null;
		}
	}

	/**
	 * Should always be synchronised on seriesCollection when called.
	 */
	public XYSeries getXYSeries(SeriesKey key)
	{
		synchronized (seriesCollection) {
			if (keyToSeries.containsKey(key)) {
				return keyToSeries.get(key);
			}

			return null;
		}
	}

	/**
	 * Should always be synchronised on seriesCollection when called.
	 * @return >0 when series found.
	 */
	public int getJFreeChartIndex(SeriesKey key)
	{
		synchronized (seriesCollection) {
			XYSeries series = keyToSeries.get(key);

			for (int i = 0; i < seriesCollection.getSeriesCount(); i++) {
				if (seriesCollection.getSeries(i).equals((series)))
					return i;
			}

			return -1;
		}
	}

	/**
	 * Getter for property graphTitle.
	 * @return Value of property graphTitle.
	 */
	public String getTitle()
	{
		return graphTitle.getStringValue();
	}

	/**
	 * Setter for property graphTitle.
	 * @param value Value of property graphTitle.
	 */
	public void setTitle(String value)
	{
		try {
			graphTitle.setValue(value);
			doEnables();
			updateGraph();
		} catch (SettingException e) {
			// Shouldn't happen.
		}
	}

	/**
	 * Getter for property titleFont.
	 * @return Value of property titleFont.
	 */
	public FontColorPair getTitleFont()
	{
		return titleFont.getFontColorValue();
	}

	/**
	 * Setter for property titleFont.
	 * @param font Value of property titleFont.
	 */
	public void setTitleFont(FontColorPair font)
	{
		try {
			titleFont.setValue(font);
			doEnables();
			updateGraph();
		} catch (SettingException e) {
			// Shouldn't happen.
		}
	}

	/**
	 * Getter for property legendFont.
	 * @return Value of property legendFont.
	 */
	public FontColorPair getLegendFont()
	{
		return legendFont.getFontColorValue();
	}

	/**
	 * Setter for property legendFont.
	 * @param font Value of property legendFont.
	 */
	public void setLegendFont(FontColorPair font)
	{
		try {
			legendFont.setValue(font);
			doEnables();
			updateGraph();
		} catch (SettingException e) {
			// Shouldn't happen.
		}
	}

	/**
	 * Getter for property legendVisible.
	 * @return Value of property legendVisible.
	 */
	public boolean isLegendVisible()
	{
		return legendVisible.getBooleanValue();
	}

	/**
	 * Setter for property legendVisible.
	 * @param visible Value of property legendVisible.
	 */
	public void setLegendVisible(boolean visible)
	{
		try {
			legendVisible.setValue(visible);
			doEnables();
			updateGraph();
		} catch (SettingException e) {
			// Shouldn't happen.
		}
	}

	/**
	 * Getter for property logarithmic.
	 * @return the legend's position index:
	 * <ul>
	 *	<li>0: LEFT
	 *	<li>1: RIGHT
	 *	<li>2: BOTTOM
	 *  <li>3: TOP
	 * </ul>
	 */
	public int getLegendPosition()
	{
		return legendPosition.getCurrentIndex();
	}

	/**
	 * Setter for property logarithmic.
	 * @param value Represents legend position
	 * <ul>
	 *	<li>0: LEFT
	 *	<li>1: RIGHT
	 *	<li>2: BOTTOM
	 *	<li>4: TOP
	 * </ul>
	 */
	public void setLegendPosition(int value) throws SettingException
	{
		legendPosition.setSelectedIndex(value);
		doEnables();
		updateGraph();
	}

	/**
	 * Return settings of the x-Axis.
	 * 
	 * @return Settings of the x-Axis.
	 */
	public AxisSettings getXAxisSettings()
	{
		return xAxisSettings;
	}

	/**
	 * Return settings of the y-Axis.
	 * 
	 * @return Settings of the y-Axis.
	 */
	public AxisSettings getYAxisSettings()
	{
		return yAxisSettings;
	}

	/**
	 * Return display settings of the graph.
	 * 
	 * @return Display settings of the graph.
	 */
	public DisplaySettings getDisplaySettings()
	{
		return displaySettings;
	}

	private String getUniqueSeriesName(String seriesName)
	{
		synchronized (seriesCollection) {
			int counter = 0;
			String name = seriesName;

			/* Name sure seriesName is unique */
			while (true) {
				boolean nameExists = false;

				for (Map.Entry<SeriesKey, XYSeries> entry : keyToSeries.entrySet()) {
					if (name.equals(entry.getValue().getKey())) {
						nameExists = true;
						break;
					}
				}

				if (nameExists) {
					counter++;
					name = seriesName + " (" + counter + ")";
				} else {
					break;
				}
			}

			return name;
		}
	}

	public void moveUp(java.util.Vector<SeriesKey> keys)
	{
		synchronized (seriesCollection) {
			XYSeries[] newOrder = new XYSeries[seriesCollection.getSeriesCount()];
			java.util.Vector<XYSeries> moveUpSet = new java.util.Vector<XYSeries>();

			for (int i = 0; i < newOrder.length; i++)
				newOrder[i] = seriesCollection.getSeries(i);

			for (SeriesKey key : keys) {
				if (keyToSeries.containsKey(key))
					moveUpSet.add(keyToSeries.get(key));
			}

			for (int i = 1; i < newOrder.length; i++) {
				if (moveUpSet.contains(newOrder[i])) {
					XYSeries tmp = newOrder[i];
					newOrder[i] = newOrder[i - 1];
					newOrder[i - 1] = tmp;
				}
			}

			XYSeriesCollection newCollection = new XYSeriesCollection();

			for (int i = 0; i < newOrder.length; i++)
				newCollection.addSeries(newOrder[i]);

			plot.setDataset(newCollection);

			this.seriesCollection = newCollection;
			this.seriesList.updateSeriesList();
		}
	}

	public void moveDown(java.util.Vector<SeriesKey> keys)
	{
		synchronized (seriesCollection) {
			XYSeries[] newOrder = new XYSeries[seriesCollection.getSeriesCount()];
			java.util.Vector<XYSeries> moveDownSet = new java.util.Vector<XYSeries>();

			for (int i = 0; i < newOrder.length; i++)
				newOrder[i] = seriesCollection.getSeries(i);

			for (SeriesKey key : keys) {
				if (keyToSeries.containsKey(key))
					moveDownSet.add(keyToSeries.get(key));
			}

			for (int i = newOrder.length - 2; i >= 0; i--) {
				if (moveDownSet.contains(newOrder[i])) {
					XYSeries tmp = newOrder[i];
					newOrder[i] = newOrder[i + 1];
					newOrder[i + 1] = tmp;
				}
			}

			XYSeriesCollection newCollection = new XYSeriesCollection();

			for (int i = 0; i < newOrder.length; i++)
				newCollection.addSeries(newOrder[i]);

			plot.setDataset(newCollection);

			this.seriesCollection = newCollection;
			this.seriesList.updateSeriesList();
		}
	}

	/**
	 * Add a series to the buffered graph data.
	 * 
	 * @param seriesName
	 *            Name of series to add to graph.
	 */
	private SeriesKey addAdditionalPointListSeries(PointList pointList)
	{
		SeriesKey key;

		synchronized (seriesCollection) {
			String seriesName = getUniqueSeriesName(pointList.getName());

			// create a new XYSeries without sorting, allowing duplicates only if multipleXallowed == true
			PrismXYSeries newSeries = new PrismXYSeries(seriesName, true);
			newSeries.setMultipleXAllowed(true);
			this.seriesCollection.addSeries(newSeries);
			// allocate a new cache for this series

			key = new SeriesKey();

			this.keyToSeries.put(key, newSeries);
			this.keyToPointList.put(key, pointList);
			this.graphCache.put(key, new LinkedList<XYDataItem>());
			this.polyCache.put(key, new LinkedList<double[]>());

			SeriesSettings graphSeries = new SeriesSettings(this, key);
			this.keyToGraphSeries.put(key, graphSeries);
			graphSeries.addObserver(this);

			this.seriesList.updateSeriesList();
		}

		return key;
	}

	/**
	 * Add a series to the buffered graph data.
	 * 
	 * @param seriesName
	 *            Name of series to add to graph.
	 */
	public SeriesKey addSeries(String seriesName)
	{

		SeriesKey key;

		synchronized (seriesCollection) {
			seriesName = getUniqueSeriesName(seriesName);

			// create a new XYSeries without sorting, allowing duplicates only if multipleXallowed == true
			PrismXYSeries newSeries = new PrismXYSeries(seriesName, false);
			newSeries.setMultipleXAllowed(false);
			this.seriesCollection.addSeries(newSeries);
			// allocate a new cache for this series

			key = new SeriesKey();

			this.keyToSeries.put(key, newSeries);
			this.graphCache.put(key, new LinkedList<XYDataItem>());
			this.polyCache.put(key, new LinkedList<double[]>());

			SeriesSettings graphSeries = new SeriesSettings(this, key);
			this.keyToGraphSeries.put(key, graphSeries);
			graphSeries.addObserver(this);

			this.seriesList.updateSeriesList();
		}

		return key;
	}

	/**
	 * Changes the name of a series.
	 * 
	 * @param key The key identifying the series.
	 * @param seriesName New name of series.
	 */
	public void changeSeriesName(SeriesKey key, String seriesName)
	{
		synchronized (seriesCollection) {
			seriesName = getUniqueSeriesName(seriesName);

			if (keyToSeries.containsKey(key)) {
				XYSeries series = keyToSeries.get(key);
				series.setKey(seriesName);
			}
		}
	}

	/**
	 * Wholly remove a series from the current graph, by key.
	 * @param seriesKey SeriesKey of series to remove.
	 */
	public void removeSeries(SeriesKey seriesKey)
	{
		synchronized (seriesCollection) {
			// Delete from keyToSeries and seriesCollection.
			if (keyToSeries.containsKey(seriesKey)) {
				XYSeries series = keyToSeries.get(seriesKey);
				seriesCollection.removeSeries(series);
				keyToSeries.remove(seriesKey);
			}

			// Remove any cache.
			if (graphCache.containsKey(seriesKey)) {
				graphCache.remove(seriesKey);
			}

			if (polyCache.containsKey(seriesKey)) {
				polyCache.remove(seriesKey);
			}

			if (keyToGraphSeries.containsKey(seriesKey)) {
				keyToGraphSeries.get(seriesKey).deleteObservers();
				keyToGraphSeries.remove(seriesKey);
			}

			this.seriesList.updateSeriesList();
		}

		seriesList.updateSeriesList();
	}

	public void addPolyToSeries(SeriesKey seriesKey, double[] poly)
	{
		synchronized (seriesCollection) {
			if (polyCache.containsKey(seriesKey)) {
				LinkedList<double[]> seriesCache = polyCache.get(seriesKey);
				seriesCache.add(poly);
			}
		}

	}

	/**
	 * Add a point to the specified graph series.
	 * @param seriesKey Key of series to update.
	 * @param dataItem XYDataItem object to insert into this series.
	 */
	public void addPointToSeries(SeriesKey seriesKey, XYDataItem dataItem)
	{
		synchronized (seriesCollection) {
			if (graphCache.containsKey(seriesKey)) {

				if (true) {
					LinkedList<XYDataItem> seriesCache = graphCache.get(seriesKey);
					seriesCache.add(dataItem);
				}
			}
		}
	}

	/**
	 * Remove all points from a graph series and its cache.
	 * 
	 * @param seriesKey
	 *            Key of series to update.
	 */
	public void removeAllPoints(SeriesKey seriesKey)
	{
		synchronized (seriesCollection) {
			if (graphCache.containsKey(seriesKey)) {
				LinkedList<XYDataItem> seriesCache = graphCache.get(seriesKey);
				seriesCache.clear();
			}

			if (polyCache.containsKey(seriesKey)) {
				LinkedList<double[]> seriesCache = polyCache.get(seriesKey);
				seriesCache.clear();
			}

			if (keyToSeries.containsKey(seriesKey)) {
				XYSeries series = keyToSeries.get(seriesKey);
				series.clear();
			}
		}
	}

	/**
	 * TODO: Document this!
	 */
	public InputSource resolveEntity(String publicId, String systemId) throws SAXException, IOException
	{
		InputSource inputSource = null;

		// override the resolve method for the dtd
		if (systemId.endsWith("dtd")) {
			// get appropriate dtd from classpath
			InputStream inputStream = Graph.class.getClassLoader().getResourceAsStream("dtds/chartformat.dtd");
			if (inputStream != null)
				inputSource = new InputSource(inputStream);
		}
		return inputSource;
	}

	/** Refactored out from load(), parses an Axis. */
	public static void parseAxis(ValueAxis axis, String minValue, String maxValue, String majorGridInterval, String minorGridInterval)
	{
		double min, max, major, minor;

		// can't work with null axis
		if (axis == null)
			return;

		try {
			min = Double.parseDouble(minValue);
		} catch (NumberFormatException e) {
			min = 0;
		}
		try {
			max = Double.parseDouble(maxValue);
		} catch (NumberFormatException e) {
			if (min < 1) {
				max = 1.0;
			} else {
				max = min + 1;
			}
		}
		try {
			major = Double.parseDouble(majorGridInterval);
		} catch (NumberFormatException e) {
			major = max / 5;
		}
		try {
			minor = Double.parseDouble(minorGridInterval);
		} catch (NumberFormatException e) {
			minor = major / 10;
		}

		// set parameters for axis
		axis.setLowerBound(min);
		axis.setUpperBound(max);
		axis.setTickMarkInsideLength((float) minor);
		axis.setTickMarkInsideLength((float) major);
	}

	public static boolean parseBoolean(String boolStr)
	{
		return ("true").equals(boolStr);
	}

	public static int parseInt(String intStr)
	{
		try {
			int d = Integer.parseInt(intStr);
			return d;
		} catch (NumberFormatException e) {
			return 0;
		}
	}

	public static double parseDouble(String doubleStr)
	{
		try {
			double d = Double.parseDouble(doubleStr);
			return d;
		} catch (NumberFormatException e) {
			return Double.NaN;
		}
	}

	/** Refactored out from load(), parses a Font */
	public static Font parseFont(String fontName, String fontStyle, String fontSize)
	{
		int style, size;
		try {
			size = Integer.parseInt(fontSize);
			style = Integer.parseInt(fontStyle);
		} catch (NumberFormatException e) {
			// If there's an error, set defaults
			size = 14;
			style = Font.PLAIN;
		}
		if (size <= 0)
			size = 12;
		if (fontName.equals(""))
			fontName = "SansSerif";

		return new Font(fontName, style, size);
	}

	/** Refactored out from load(), parses a Color. */
	public static Color parseColor(String red, String green, String blue)
	{
		int r, g, b;

		try {
			r = Integer.parseInt(red);
			g = Integer.parseInt(green);
			b = Integer.parseInt(blue);

			if (r > 255)
				r = 255;
			if (r < 0)
				r = 0;
			if (g > 255)
				g = 255;
			if (g < 0)
				g = 0;
			if (b > 255)
				b = 255;
			if (b < 0)
				b = 0;
		} catch (NumberFormatException e) {
			// If theres an error, set defaults
			r = 0;
			g = 0;
			b = 0;
		}
		return new Color(r, g, b);
	}

	/**
	 * Method to load a PRISM 'gra' file into the application.
	 * @param file Name of the file to load.
	 * @return The model of the graph contained in the file.
	 * @throws GraphException if I/O errors have occurred.
	 */
	public static Graph load(File file) throws GraphException
	{

		Graph graph = new Graph();

		try {
			DocumentBuilderFactory factory = DocumentBuilderFactory.newInstance();

			factory.setValidating(true);
			factory.setIgnoringElementContentWhitespace(true);
			DocumentBuilder builder = factory.newDocumentBuilder();
			builder.setEntityResolver(graph);
			Document doc = builder.parse(file);
			Element chartFormat = doc.getDocumentElement();

			graph.setTitle(chartFormat.getAttribute("graphTitle"));

			String titleFontName = chartFormat.getAttribute("titleFontName");
			String titleFontSize = chartFormat.getAttribute("titleFontSize");
			String titleFontStyle = chartFormat.getAttribute("titleFontStyle");

			Font titleFont = parseFont(titleFontName, titleFontStyle, titleFontSize);

			String titleFontColourR = chartFormat.getAttribute("titleFontColourR");
			String titleFontColourG = chartFormat.getAttribute("titleFontColourG");
			String titleFontColourB = chartFormat.getAttribute("titleFontColourB");
			Color titleFontColour = parseColor(titleFontColourR, titleFontColourG, titleFontColourB);

			graph.setTitleFont(new FontColorPair(titleFont, titleFontColour));
			graph.setLegendVisible(parseBoolean(chartFormat.getAttribute("legendVisible")));

			String legendPosition = chartFormat.getAttribute("legendPosition");

			// Facilitate for bugs export in previous prism versions.
			if (chartFormat.getAttribute("versionString").equals(""))
				graph.setLegendPosition(RIGHT);
			else {
				if (legendPosition.equals("left"))
					graph.setLegendPosition(LEFT);
				else if (legendPosition.equals("right"))
					graph.setLegendPosition(RIGHT);
				else if (legendPosition.equals("bottom"))
					graph.setLegendPosition(BOTTOM);
				else if (legendPosition.equals("top"))
					graph.setLegendPosition(TOP);
				else
					// Probably was manual, now depricated
					graph.setLegendPosition(RIGHT);
			}

			//Get the nodes used to describe the various parts of the graph
			NodeList rootChildren = chartFormat.getChildNodes();

			// Element layout is depricated for now. 
			Element layout = (Element) rootChildren.item(0);
			Element xAxis = (Element) rootChildren.item(1);
			Element yAxis = (Element) rootChildren.item(2);

			graph.getXAxisSettings().load(xAxis);
			graph.getYAxisSettings().load(yAxis);

			//Read the headings and widths for each series 
			for (int i = 3; i < rootChildren.getLength(); i++) {
				Element series = (Element) rootChildren.item(i);
				SeriesKey key = graph.addSeries(series.getAttribute("seriesHeading"));

				synchronized (graph.getSeriesLock()) {
					SeriesSettings seriesSettings = graph.getGraphSeries(key);
					seriesSettings.load(series);

					NodeList graphChildren = series.getChildNodes();

					//Read each series out of the file and add its points to the graph
					for (int j = 0; j < graphChildren.getLength(); j++) {
						Element point = (Element) graphChildren.item(j);
						graph.addPointToSeries(key, new XYDataItem(parseDouble(point.getAttribute("x")), parseDouble(point.getAttribute("y"))));
					}
				}
			}

			//Return the model of the graph 
			return graph;
		} catch (Exception e) {
			throw new GraphException("Error in loading chart: " + e);
		}
	}

	/**
	 * Allows graphs to be saved to the PRISM 'gra' file format.
	 * 
	 * @param file
	 *            The file to save the graph to.
	 * @throws GraphException
	 *             If the file cannot be written.
	 */
	public void save(File file) throws PrismException
	{
		try {
			JFreeChart chart = getChart();
			DocumentBuilderFactory factory = DocumentBuilderFactory.newInstance();
			DocumentBuilder builder = factory.newDocumentBuilder();
			Document doc = builder.newDocument();

			Element chartFormat = doc.createElement("chartFormat");

			chartFormat.setAttribute("versionString", prism.Prism.getVersion());
			chartFormat.setAttribute("graphTitle", getTitle());

			Font titleFont = getTitleFont().f;
			chartFormat.setAttribute("titleFontName", titleFont.getName());
			chartFormat.setAttribute("titleFontSize", "" + titleFont.getSize());
			chartFormat.setAttribute("titleFontStyle", "" + titleFont.getStyle());

			Color titleFontColor = (Color) getTitleFont().c;
			chartFormat.setAttribute("titleFontColourR", "" + titleFontColor.getRed());
			chartFormat.setAttribute("titleFontColourG", "" + titleFontColor.getGreen());
			chartFormat.setAttribute("titleFontColourB", "" + titleFontColor.getBlue());

			chartFormat.setAttribute("legendVisible", isLegendVisible() ? "true" : "false");

			Font legendFont = getLegendFont().f;
			chartFormat.setAttribute("legendFontName", "" + legendFont.getName());
			chartFormat.setAttribute("legendFontSize", "" + legendFont.getSize());
			chartFormat.setAttribute("legendFontStyle", "" + legendFont.getStyle());

			Color legendFontColor = getLegendFont().c;

			chartFormat.setAttribute("legendFontColourR", "" + legendFontColor.getRed());
			chartFormat.setAttribute("legendFontColourG", "" + legendFontColor.getGreen());
			chartFormat.setAttribute("legendFontColourB", "" + legendFontColor.getBlue());

			switch (getLegendPosition()) {
			case LEFT:
				chartFormat.setAttribute("legendPosition", "left");
				break;
			case BOTTOM:
				chartFormat.setAttribute("legendPosition", "bottom");
				break;
			case TOP:
				chartFormat.setAttribute("legendPosition", "top");
				break;
			default:
				chartFormat.setAttribute("legendPosition", "right");
			}

			Element layout = doc.createElement("layout");
			chartFormat.appendChild(layout);

			Element xAxis = doc.createElement("axis");
			getXAxisSettings().save(xAxis);
			chartFormat.appendChild(xAxis);

			Element yAxis = doc.createElement("axis");
			getYAxisSettings().save(yAxis);
			chartFormat.appendChild(yAxis);

			synchronized (getSeriesLock()) {
				/* Make sure we preserve ordering. */
				for (int i = 0; i < seriesList.getSize(); i++) {
					SeriesKey key = seriesList.getKeyAt(i);

					Element series = doc.createElement("graph");
					SeriesSettings seriesSettings = getGraphSeries(key);
					seriesSettings.save(series);

					XYSeries seriesData = getXYSeries(key);

					for (int j = 0; j < seriesData.getItemCount(); j++) {
						Element point = doc.createElement("point");

						point.setAttribute("x", "" + seriesData.getX(j));
						point.setAttribute("y", "" + seriesData.getY(j));

						series.appendChild(point);
					}

					chartFormat.appendChild(series);
				}
			}

			doc.appendChild(chartFormat);

			// File writing 
			Transformer t = TransformerFactory.newInstance().newTransformer();
			t.setOutputProperty("doctype-system", "chartformat.dtd");
			t.setOutputProperty("indent", "yes");
			t.transform(new DOMSource(doc), new StreamResult(new FileOutputStream(file)));
		} catch (IOException e) {
			throw new PrismException(e.getMessage());
		} catch (DOMException e) {
			throw new PrismException("Problem saving graph: DOM Exception: " + e);
		} catch (ParserConfigurationException e) {
			throw new PrismException("Problem saving graph: Parser Exception: " + e);
		} catch (TransformerConfigurationException e) {
			throw new PrismException("Problem saving graph: Error in creating XML: " + e);
		} catch (TransformerException e) {
			throw new PrismException("Problem saving graph: Transformer Exception: " + e);
		} catch (SettingException e) {
			throw new PrismException(e.getMessage());
		}
	}

	/**
	 * Exports the current graph to Matlab file format.
	 * @param f The file to write the data to.
	 */
	public void exportToMatlab(File f) throws IOException
	{
		PrintWriter out = new PrintWriter(new FileWriter(f));

		out.println("%=========================================");
		out.println("%Generated by PRISM Chart Package");
		out.println("%=========================================");
		out.println();

		//Seriesdata
		synchronized (getSeriesLock()) {
			/* Make sure we preserve ordering. */
			for (int i = 0; i < seriesList.getSize(); i++) {
				StringBuffer x = new StringBuffer("x" + i + " = [");
				StringBuffer y = new StringBuffer("y" + i + " = [");

				SeriesKey key = seriesList.getKeyAt(i);

				XYSeries seriesData = getXYSeries(key);

				for (int j = 0; j < seriesData.getItemCount(); j++) {
					x.append(seriesData.getX(j) + " ");
					y.append(seriesData.getY(j) + " ");
				}

				x.append("];");
				y.append("];");

				out.println(x.toString());
				out.println(y.toString());
			}

			//Create a figure
			out.println();
			out.println("figure1 = figure('Color', [1 1 1], 'PaperPosition',[0.6345 6.345 20.3 15.23],'PaperSize',[20.98 29.68]);");

			//Create axes
			boolean xLog = getXAxisSettings().isLogarithmic();
			boolean yLog = getYAxisSettings().isLogarithmic();

			out.println();

			if (xLog && yLog)
				out.println("axes1 = axes('Parent', figure1, 'FontSize', 16, 'XScale', 'log', 'YScale', 'log');");
			else if (xLog)
				out.println("axes1 = axes('Parent', figure1, 'FontSize', 16, 'XScale', 'log');");
			else if (yLog)
				out.println("axes1 = axes('Parent', figure1, 'FontSize', 16, 'YScale', 'log');");
			else
				out.println("axes1 = axes('Parent', figure1, 'FontSize', 16);");

			out.println("xlabel(axes1, '" + getXAxisSettings().getHeading() + "');");
			out.println("ylabel(axes1, '" + getYAxisSettings().getHeading() + "');");

			out.println("box(axes1, 'on');");
			out.println("hold(axes1, 'all');");

			//Graph title
			out.println();

			String title = "";
			StringTokenizer st = new StringTokenizer(getTitle(), "\n");

			int num = st.countTokens();
			for (int i = 0; i < num; i++) {
				title += "'" + st.nextToken() + "'";
				if (i < num - 1)
					title += ", char(10),";
			}

			out.println("title(axes1,[" + title + "])");

			//Sort out logarithmic scales
			String scaleType = "plot";
			if (seriesList.getSize() > 0) {
				if (xLog && yLog)
					scaleType = "loglog";
				else if (xLog)
					scaleType = "semilogx";
				else if (yLog)
					scaleType = "semilogy";
			}

			//Create plots
			for (int i = 0; i < seriesList.getSize(); i++) {
				SeriesKey key = seriesList.getKeyAt(i);
				SeriesSettings seriesSettings = getGraphSeries(key);

				String marker = "'";
				if (seriesSettings.showPoints() && seriesSettings.getSeriesShape() != SeriesSettings.NONE) {
					switch (seriesSettings.getSeriesShape()) {
					case SeriesSettings.CIRCLE:
						marker += "o";
						break;
					case SeriesSettings.SQUARE:
						marker += "s";
						break;
					case SeriesSettings.TRIANGLE:
						marker += "^";
						break;
					case SeriesSettings.RECTANGLE_H:
						marker += "d";
						break;
					case SeriesSettings.RECTANGLE_V:
						marker += "x";
						break;
					}
				}

				if (seriesSettings.showLines()) {
					switch (seriesSettings.getLineStyle()) {
					case SeriesSettings.SOLID:
						marker += "-";
						break;
					case SeriesSettings.DASHED:
						marker += "--";
						break;
					case SeriesSettings.DOT_DASHED:
						marker += "-.";
						break;
					}
				}
				marker += "'";
				out.println("plot" + i + " = " + scaleType + "(x" + i + ", y" + i + ", " + marker + ", 'Parent', axes1, 'LineWidth', 2);");
			}

			//			Create legend
			String seriesNames = "";
			for (int i = 0; i < seriesList.getSize(); i++) {
				SeriesKey key = seriesList.getKeyAt(i);
				SeriesSettings seriesSettings = getGraphSeries(key);

				seriesNames += "'" + seriesSettings.getSeriesHeading() + "'";
				if (i < seriesList.getSize() - 1)
					seriesNames += ", ";
			}

			//Determine location

			String loc = "";
			switch (legendPosition.getCurrentIndex()) {
			case LEFT:
				loc = "'WestOutside'";
				break;
			case RIGHT:
				loc = "'EastOutside'";
				break;
			case BOTTOM:
				loc = "'SouthOutside'";
				break;
			case TOP:
				loc = "'NorthOutside'";
				break;
			}

			if (isLegendVisible())
				out.println("legend1 = legend(axes1,{" + seriesNames + "},'Location', " + loc + ");");

			out.flush();
			out.close();
		}
	}

	/**
	 * Renders the current graph to a JPEG file.
	 * 
	 * @param file
	 *             The file to export the JPEG data to.
	 * @throws GraphException, IOException
	 *             If file cannot be written to.
	 */
	public void exportToJPEG(File file, int width, int height) throws GraphException, IOException
	{
		ChartUtilities.saveChartAsJPEG(file, 1.0f, this.chart, width, height);
	}

	public void exportToEPS(File file, int width, int height) throws GraphException, IOException
	{
		FileOutputStream fileOutputStream = new FileOutputStream(file);

		EpsGraphics g2d = new EpsGraphics(this.getTitle(), fileOutputStream, 0, 0, width, height, ColorMode.COLOR_RGB);

		// Don't export fonts as vectors, no hope of getting same font as publication.
		g2d.setAccurateTextMode(false); // Does not rotate y-axis label.

		chart.draw(g2d, new Rectangle(width, height));

		g2d.close();
		g2d.dispose();
	}

	/**
	 * Renders the current graph to a JPEG file.
	 * 
	 * @param file
	 *             The file to export the JPEG data to.
	 * @throws GraphException, IOException
	 *             If file cannot be written to.
	 */
	public void exportToPNG(File file, int width, int height, boolean alpha) throws GraphException, IOException
	{

		FileOutputStream fileOutputStream = new FileOutputStream(file);

		KeypointPNGEncoderAdapter encoder = new KeypointPNGEncoderAdapter();
		encoder.setEncodingAlpha(alpha);

		Paint bgPaint = chart.getBackgroundPaint();

		if (alpha) {
			chart.setBackgroundPaint(null);
		}

		BufferedImage bufferedImage = chart.createBufferedImage(width, height, alpha ? BufferedImage.TYPE_INT_ARGB : BufferedImage.TYPE_INT_RGB, null);

		if (alpha) {
			chart.setBackgroundPaint(bgPaint);
		}

		encoder.encode(bufferedImage, fileOutputStream);

		fileOutputStream.flush();
		fileOutputStream.close();

		//ChartUtilities.saveChartAsPNG(file, this.chart, width, height, null, alpha, 9);
	}

	/**
	 * Exports the current data sets to a comma seperated value file.
	 * 
	 * @param file
	 *            The file to export the csv data to.
	 * @throws GraphException
	 *             If file cannont be written to.
	 */
	public void exportToCSV(File file) throws GraphException
	{
		/*try 
		{ 
			PrintWriter out = new PrintWriter(new FileWriter(file)); 
			
			// Create the header, this is a comma seperate list of the titles 
			// of the series followed by an x or y String header = ""; 
			  
			synchronize(this.getSeriesLock())
			{
				
			}
			
			int largestGraph = 0;
			 
			// Consider every series in this model 
			for (int i = 0; i < getNumSeries(); i++) 
			{ 
				// While were iterating through the graphs, record the size of the largest graph 
				if (getGraphPoints(i).getItemCount() > largestGraph) 
				{ 
					largestGraph = getGraphPoints(i).getItemCount(); 
				}
			   
				// Add the series titles 
				header += "[" + getGraphPoints(i).getKey() + "].x,"; 
				header += "[" + getGraphPoints(i).getKey() + "].y";
			 
				// If this isnt the last series add the comma 
				if (i < getNumSeries() - * 1) 
					header += ", "; 
			}
			
			// Print the header to the file
		    out.println(header);
		 
		    // Format of the file has the data in columns, so iterate through 
		    // the points in all the graphs adding the next point from each graph 
		    // on the same line 
		    
		    String line = ""; 
		  
		    for (int i = 0; i < largestGraph; i++) 
		    { 
		    	line = "";
		  
		    	// For each point, iterate through all the graphs 
		    	for (int j = 0; j < getNumSeries(); j++) 
		    	{ 
		    		XYSeries gs = getGraphPoints(j);
		    		  
		    		// If the current graph doesnt have this point, ie its smaller 
		    		// than the largest supply "","" for the data 
		  
					if (i >= gs.getItemCount()) 
					{ 
						line += "\"\",\"\""; 
					} 
					else 
					{ 
						// Otherwise add the next data point 
						line += gs.getX(i) + ","; line += gs.getY(i); 
					}
		  
					if (j < getNumSeries() - 1) 
						line += ", "; 
		    	}
		 
		    	out.println(line); 
		    }
		    
		    out.flush(); 
		    out.close();
		} 
		catch (Exception e) 
		{ 
			throw new ChartException(e); 
		}	
		*/
	}

	/**
	 * This inner class provides a means of asynchronously performing graph
	 * updates. This has two advantages.
	 * 
	 * Firstly, it improves the callee's performance as adding data is
	 * effectively non-blocking.
	 * 
	 * Secondly, it improves the interactive response of the chart, as fewer
	 * draws are performed. According to JFreeChart's own FAQ, it is not
	 * designed for real-time charting.
	 */
	private class GraphUpdateTask extends TimerTask
	{
		private void processGraphCache(HashMap<SeriesKey, LinkedList<XYDataItem>> graphCache)
		{
			synchronized (seriesCollection) {
				for (Map.Entry<SeriesKey, LinkedList<XYDataItem>> entry : graphCache.entrySet()) {
					/* The series key should map to a series. */
					if (keyToSeries.containsKey(entry.getKey())) {
						XYSeries series = keyToSeries.get(entry.getKey());
						Color c = keyToGraphSeries.get(entry.getKey()).getSeriesColour();
						LinkedList<XYDataItem> seriesCache = entry.getValue();
						while (!seriesCache.isEmpty()) {
							XYDataItem item = seriesCache.removeFirst();
							series.addOrUpdate(item.getX(), item.getY());
						}
					}
				}
				// draw polys as annotations
				for (Map.Entry<SeriesKey, LinkedList<double[]>> entry : polyCache.entrySet()) {

					/* The series key should map to a series. */
					if (keyToSeries.containsKey(entry.getKey())) {
						LinkedList<double[]> seriesCache = entry.getValue();
						while (!seriesCache.isEmpty()) {
							double[] poly = seriesCache.removeFirst();
							Color c = keyToGraphSeries.get(entry.getKey()).getSeriesColour();
							Color c_new = new Color(c.getRed(), c.getGreen(), c.getBlue(), 150); // make translucent
							XYPolygonAnnotation ann = new XYPolygonAnnotation(poly, new BasicStroke(), Color.BLACK, c_new);
							plot.addAnnotation(ann);
						}
					}
				}

			}
		}

		public void run()
		{
			processGraphCache(graphCache);
		}
	}

	/**
	 * This is a dummy class used to be able to index Series in this graph.
	 * Previously this was done using integers, which was unsafe if removeSeries
	 * was used. The hashcode() and equals() implementation of Object (based on
	 * object identity) are sufficient to use this as the key of a HashMap.
	 * In addition, we add a 'next' field, to allow the same class to be used
	 * to store (linked) lists of keys.
	 */
	public class SeriesKey
	{
		public SeriesKey next = null;
	}

}
