//==============================================================================
//
//	Copyright (c) 2026-
//	Authors:
//	* Gabriel Santos <gabriel.santos@cs.ox.ac.uk> (University of Oxford)
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

package explicit;

import java.util.BitSet;

/**
 * The objectives of the coalitions in multi-player equilibria (more than two coalitions), and when each objective is
 * decided: a coalition is done (D) once its target is reached (or, for rewards, when it has nothing left to collect),
 * failed (E) once its objective can no longer be satisfied (until violated, or bound expired), and undecided otherwise.
 * Shared by the model checking ({@link CSGModelCheckerEquilibria}) and the strategies it generates (to update their
 * memory), so that both follow exactly the same rules.
 */
public class CSGMultiObjectives
{
	/** Kinds of objectives */
	public static final int OBJ_UNBOUNDED = 0, OBJ_BOUNDED = 1, OBJ_NEXT = 2, OBJ_CUMULATIVE = 3, OBJ_INSTANTANEOUS = 4;

	/** Number of coalitions */
	public final int n;
	/** Whether the objectives are rewards (otherwise probabilities) */
	public final boolean rew;
	/** Kind (OBJ_*) and bound (-1 if unbounded) of each coalition's objective */
	public final int[] kind, bound;
	/** Target states (null for cumulative/instantaneous rewards) and states to remain in (until, or null) */
	public final BitSet[] targets, remain;
	/** Coalitions with bounded objectives, and the largest bound (the horizon of the bounded phase) */
	public final BitSet bounded = new BitSet();
	public final int horizon;

	public CSGMultiObjectives(int n, boolean rew, int[] kind, int[] bound, BitSet[] targets, BitSet[] remain)
	{
		this.n = n;
		this.rew = rew;
		this.kind = kind;
		this.bound = bound;
		this.targets = targets;
		this.remain = remain;
		int h = 0;
		for (int c = 0; c < n; c++) {
			if (kind[c] != OBJ_UNBOUNDED) {
				bounded.set(c);
				h = Math.max(h, bound[c]);
			}
		}
		horizon = h;
	}

	/**
	 * Status of (not yet decided) coalition c in state s after t steps: 1 done, -1 failed, 0 undecided.
	 * t is only relevant to bounded objectives (which are all decided once t reaches the horizon).
	 */
	public int status(int c, int s, int t)
	{
		boolean target = targets[c] != null && targets[c].get(s);
		boolean inRemain = remain == null || remain[c] == null || remain[c].get(s);
		switch (kind[c]) {
		case OBJ_UNBOUNDED:
			return target ? 1 : (!rew && !inRemain) ? -1 : 0;
		case OBJ_BOUNDED:
			return (target && t <= bound[c]) ? 1 : (!inRemain || t >= bound[c]) ? -1 : 0;
		case OBJ_NEXT:
			return t == 0 ? 0 : target ? 1 : -1;
		case OBJ_CUMULATIVE:
		case OBJ_INSTANTANEOUS:
			return t >= bound[c] ? 1 : 0;
		default:
			return 0;
		}
	}

	/**
	 * Adds to D and E the coalitions (not already in either) that become done or fail in state s after t steps.
	 */
	public void classify(int s, int t, BitSet D, BitSet E)
	{
		BitSet newD = new BitSet(), newE = new BitSet();
		for (int c = 0; c < n; c++) {
			if (D.get(c) || E.get(c))
				continue;
			int st = status(c, s, t);
			if (st > 0)
				newD.set(c);
			else if (st < 0)
				newE.set(c);
		}
		D.or(newD);
		E.or(newE);
	}

	/**
	 * Whether all bounded objectives are decided (the coalitions with bounded objectives are all in D or E):
	 * from then on, the step count no longer matters.
	 */
	public boolean boundedDecided(BitSet D, BitSet E)
	{
		BitSet undecided = (BitSet) bounded.clone();
		undecided.andNot(D);
		undecided.andNot(E);
		return undecided.isEmpty();
	}
}
