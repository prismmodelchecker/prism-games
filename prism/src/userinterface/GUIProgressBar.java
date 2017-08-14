//==============================================================================
//	
//	Copyright (c) 2002-
//	Authors:
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

package userinterface;

import javax.swing.*;
import java.awt.event.ActionListener;
import java.awt.event.ActionEvent;
import java.awt.Insets;
import prism.Prism;


/** This class is really for convienience.  It extends JPanel.  It sets up desirable properties for a
 * progress bar at the bottom of the screen and also provides a button to cancel an operation.
 */
public class GUIProgressBar extends JPanel
{
    JProgressBar progress;
    JButton cancel_computation;

    //CONSTANTS
    public static final int MINIMUM_WIDTH = 45;
    public static final int MINIMUM_HEIGHT = 20;
    //Constructor
    
    /** Creates a new GUIProgressBar. */
    public GUIProgressBar(int min, int max, final GUIPrism prism)
    {
	super();

	progress = new JProgressBar(min, max);
	progress.setBorder(null);
	//Progress.setHorizontalAlignment(javax.swing.SwingConstants.LEFT);
	//setBorder(new javax.swing.border.EtchedBorder());
	progress.setMinimumSize(new java.awt.Dimension(MINIMUM_WIDTH, MINIMUM_HEIGHT));
	this.add(progress);

	cancel_computation = new JButton();
	cancel_computation.setToolTipText("Cancel current computation.");
	cancel_computation.setMinimumSize(new java.awt.Dimension(MINIMUM_HEIGHT, MINIMUM_HEIGHT));
	cancel_computation.setIcon(GUIPrism.getIconFromImage("tinyStop.png"));
	cancel_computation.setMargin(new Insets(0, 0, 0, 0));
	cancel_computation.addActionListener(new ActionListener() {
		@Override
		public void actionPerformed(ActionEvent e) {
		    prism.getPrism().setCancel(true);
		}
	    });
	this.add(cancel_computation);
	cancel_computation.setEnabled(false);
    }

    public void setIndeterminate(boolean newValue)
    {
	progress.setIndeterminate(newValue);
	cancel_computation.setEnabled(newValue);
    }
    
}
