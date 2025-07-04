//==============================================================================
//	
//	Copyright (c) 2002-
//	Authors:
//	* Dave Parker <d.a.parker@cs.bham.ac.uk> (University of Birmingham/Oxford)
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

package prism;

public enum ModelType
{
	// List of model types (ordered alphabetically)
	CTMC("continuous-time Markov chain") {
		@Override
		public boolean choicesSumToOne()
		{
			return false;
		}

		@Override
		public boolean continuousTime()
		{
			return true;
		}

		@Override
		public boolean nondeterministic()
		{
			return false;
		}

		@Override
		public String probabilityOrRate()
		{
			return RATE;
		}
	},
	CTMDP("continuous-time Markov decision process") {
		@Override
		public boolean choicesSumToOne()
		{
			return false;
		}

		@Override
		public boolean continuousTime()
		{
			return true;
		}

		@Override
		public String probabilityOrRate()
		{
			return RATE;
		}

		@Override
		public ModelType removeNondeterminism()
		{
			return CTMC;
		}
	},
	/*** ***/
	CSG("concurrent stochastic game") {
		@Override
		public boolean multiplePlayers()
		{
			return true;
		}

		@Override
		public boolean concurrent()
		{
			return true;
		}

		@Override
		public ModelType removeNondeterminism()
		{
			return DTMC;
		}
	}
	,
	/*** ***/
	DTMC("discrete-time Markov chain") {
		@Override
		public boolean nondeterministic()
		{
			return false;
		}
	},
	LTS("labelled transition system") {
		@Override
		public boolean isProbabilistic()
		{
			return false;
		}

		@Override
		public String probabilityOrRate()
		{
			return NEITHER;
		}

		@Override
		public ModelType removeNondeterminism()
		{
			return DTMC;
		}
	},
	MDP("Markov decision process") {
		@Override
		public ModelType removeNondeterminism()
		{
			return DTMC;
		}
	},
	POMDP("partially observable Markov decision process") {
		@Override
		public boolean partiallyObservable()
		{
			return true;
		}
		
		@Override
		public ModelType removeNondeterminism()
		{
			return DTMC;
		}
	},
	POPTA("partially observable probabilistic timed automaton") {
		@Override
		public boolean continuousTime()
		{
			return true;
		}
		
		@Override
		public boolean realTime()
		{
			return true;
		}
		
		@Override
		public boolean partiallyObservable()
		{
			return true;
		}
		
		@Override
		public ModelType removeNondeterminism()
		{
			return DTMC;
		}
	},
	PTA("probabilistic timed automaton") {
		@Override
		public boolean continuousTime()
		{
			return true;
		}
		
		@Override
		public boolean realTime()
		{
			return true;
		}
	},
	STPG("stochastic two-player game") {
		@Override
		public boolean multiplePlayers()
		{
			return true;
		}

		@Override
		public ModelType removeNondeterminism()
		{
			return DTMC;
		}
	},
	SMG("stochastic multi-player game") {
		@Override
		public boolean multiplePlayers()
		{
			return true;
		}

		@Override
		public ModelType removeNondeterminism()
		{
			return DTMC;
		}
	},
	TPTG("turn-based probabilistic timed game") {
		@Override
		public boolean continuousTime()
		{
			return true;
		}

		@Override
		public boolean realTime()
		{
			return true;
		}

		@Override
		public boolean multiplePlayers()
		{
			return true;
		}

		@Override
		public ModelType removeNondeterminism()
		{
			return DTMC;
		}
	},
	IDTMC("interval discrete-time Markov chain") {
		@Override
		public boolean nondeterministic()
		{
			// NB: we distinguish between nondeterminism and uncertainty
			return false;
		}

		@Override
		public boolean uncertain()
		{
			return true;
		}

		@Override
		public boolean intervals()
		{
			return true;
		}
	},
	UDTMC("uncertain discrete-time Markov chain") {
		@Override
		public boolean nondeterministic()
		{
			// NB: we distinguish between nondeterminism and uncertainty
			return false;
		}

		@Override
		public boolean uncertain()
		{
			return true;
		}
	},
	IMDP("interval Markov decision process") {
		@Override
		public ModelType removeNondeterminism()
		{
			return DTMC;
		}

		@Override
		public boolean uncertain()
		{
			return true;
		}

		@Override
		public boolean intervals()
		{
			return true;
		}
	},
	UMDP("uncertain Markov decision process") {
		@Override
		public ModelType removeNondeterminism()
		{
			return DTMC;
		}

		@Override
		public boolean uncertain()
		{
			return true;
		}
	};

	private static final String PROBABILITY = "Probability";
	private static final String RATE = "Rate";
	private static final String NEITHER = "";

	private final String fullName;

	ModelType(final String fullName) {
		this.fullName = fullName;
	}

	/**
	 * Get the full name, in words, of the this model type.
	 */
	public String fullName()
	{
		return fullName;
	}

	/**
	 * Get the PRISM keyword for this model type.
	 */
	public String keyword()
	{
		return this.name().toLowerCase();
	}

	/**
	 * Do the transitions in a choice sum to 1 for this model type?
	 * Can also use this to test whether models uses rates or probabilities.
	 */
	public boolean choicesSumToOne()
	{
		return true;
	}

	/**
	 * Are time delays continuous for this model type?
	 */
	public boolean continuousTime()
	{
		return false;
	}

	/**
	 * Is this a "real time" (timed automata like) model type?
	 */
	public boolean realTime()
	{
		return false;
	}

	/**
	 * Does this model allow nondeterministic choices?
	 */
	public boolean nondeterministic()
	{
		return true;
	}

	/**
	 * Does this model have more than 1 player?
	 */
	public boolean multiplePlayers()
	{
		return false;
	}

	/**
	 * For game models, do the players make choices in a concurrent
	 * (rather than turn-based) fashion? Note that this can only be true
	 * for games so,if true, this also implies that multiplePlayers() is true.
	 */
	public boolean concurrent()
	{
		return false;
	}

	/**
	 * Is this model probabilistic?
	 */
	public boolean isProbabilistic()
	{
		return true;
	}

	/**
	 * Does this model have probabilities or rates?
	 * Returns "Probability" or "Rate", accordingly (or "" if there are neither)
	 */
	public String probabilityOrRate()
	{
		return PROBABILITY;
	}

	/**
	 * Does the model feature partial observability? 
	 */
	public boolean partiallyObservable()
	{
		return false;
	}
	
	/**
	 * Does this model have uncertainty (e.g., intervals)?
	 */
	public boolean uncertain()
	{
		return false;
	}

	/**
	 * Does this model define probabilities using intervals?
	 */
	public boolean intervals()
	{
		return false;
	}

	/**
	 * Return the model type that results from removing the nondeterminism
	 * in this model type.
	 * <br>
	 * If there is no nondeterminism (or the removal of nondeterminism is not supported),
	 * returns the same model type.
	 */
	public ModelType removeNondeterminism()
	{
		// default: same model type
		return this;
	}

	public static ModelType parseName(String name)
	{
		try {
			return valueOf(name.toUpperCase());
		} catch (IllegalArgumentException e) {
			return null;
		}
	}
}