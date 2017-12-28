// ****************************
//
// general utility functions
//
// ****************************

// libraries
import std.stdio; 			// I/O and file system
import std.math; 			// math function library
import std.algorithm; 		// algorithm for working with ranges, sequences, etc.
import std.array; 			// array operations
import std.string; 			// string function support
import std.conv; 			// automatic conversions between types
import std.typecons; 		// misc. meta-programming support, including type constructors, tuples, etc.
import std.range; 			// enables range functions
import std.mathspecial; 	// for inverse normal distribution function
import std.random; 			// for random number generation


// *** functions ***


double[] Normalize(double[] A){
	// return a normalized array corresponding to "A" 
	double[] normA;
	normA.length = A.length;
	normA[] = A[]/sum(A);
	return normA;
	}

int ArrayIndex(T)(T[] A, T x){
	// step through array A and return index corresponding to first encounter of x
	// note that arguments are templated since both boolean and int types use this function
	int index = -1; 			// default value, if there is no match
	for (int i = 0; i < A.length; ++i){
		if (A[i] == x){
			index = i;
			break;
			}
		}
	return index;
	}

int MaxIntArray(int[] A){
	// compute the maxium value of an array, regardless of type
	int maxVal = A[0];
	for (int i = 0; i < A.length; ++i){maxVal = max(maxVal, A[i]);}
	return maxVal;
	}

double[] CumSum(double[] A){
	double[] sumA;
	// calculate a running some through each element in array A
	sumA ~= A[0];
	for (int i = 1; i < A.length; ++i){sumA ~= sumA[i-1] + A[i];}
	return sumA;
	}	
	
int SelectRand(double[] A){
    // make a random selection from array A, weighted by value
	int matchSelect;
	double r, prob;
	double[] propA;
	double[] cumA;
	propA.length = A.length;
    propA[] = A[]/sum(A); 			// normalize A
    cumA = CumSum(propA); 			// cumulative sum of normalized A
    r = uniform01;
	for (int i = 0; i < cumA.length; ++i){
		if (cumA[i] > r){
			matchSelect = i;
			break;
			}
		}
    return matchSelect;
	}
	
double[] Blend(double f, double[] comp1, double[] comp2){
	// compute a weighted average blended composition from a binary mixture; "f" is the proportion of component 1
	double[] mixed;
	for (int i = 0; i < comp1.length; ++i){mixed ~= f*comp1[i] + (1-f)*comp2[i];}
	return mixed;
	}