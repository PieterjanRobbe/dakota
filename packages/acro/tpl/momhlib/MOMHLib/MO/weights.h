/*
Multiple Objective MetaHeuristics Library in C++ MOMHLib++
Copyright (C) 2001 Andrzej Jaszkiewicz.

This library is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation (www.gnu.org); 
either version 2.1 of the License, or (at your option) any later 
version.

This library is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public
License along with this library; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/

#if !defined(AFX_TWEIGHTS_H__E315A671_35A4_11D2_BA26_08002B9A2EA9__INCLUDED_)
#define AFX_TWEIGHTS_H__E315A671_35A4_11D2_BA26_08002B9A2EA9__INCLUDED_

#include "global.h"

class TPoint;

/** Class corresponding to a weight vector */
class TWeightVector  : public vector<double>
{
public:
	/** Calculates Eucildean distance between two weight vectors */
	double Distance (const TWeightVector& WeightVector);

	/** Rescales weights */
	void Rescale (TPoint& IdealPoint, TPoint& NadirPoint);

	/** Copy operator */
	TWeightVector& operator = (const TWeightVector& SourceWeights);

	/** Constructs new vector allocating space for NumberOfObjectives elements */
	TWeightVector () {
		if (NumberOfObjectives == 0) {
			cout << "TWeightVector::TWeightVector ()\n";
			cout << "NumberOfObjectives == 0\n";
			exit (1);
		}
		resize (NumberOfObjectives, 0);
	}

	/** Normalizes the weight vector such that its element sum up to 1 */
	void Normalize ();
};


/** Set of weight vectors */
class TWeightsSet : public vector<TWeightVector>
{
	/** Uniformly distributed set of weight vector generated by 
	*	GenerateUniformCover method */
	vector<TWeightVector> UniformWeightVectors;
public:

	/** Generates uniformly distributed normalized weight vectors
	*	with density defined by parameter Steps - number of different values 
	*	that may be taken by each individual weights */
	void GenerateUniformCover (int Steps);

	/** Rescales all weight vectors
	*
	*	Uses UniformWeightVectors, i.e. you may call the method
	*	many times */
	void Rescale (TPoint& IdealPoint, TPoint& NadirPoint);

	void NormalizeForEvaluation();
};

/**	Generates andom normalized wieght vector */
TWeightVector GetRandomWeightVector ();

#endif // !defined(AFX_TWEIGHTS_H__E315A671_35A4_11D2_BA26_08002B9A2EA9__INCLUDED_)
