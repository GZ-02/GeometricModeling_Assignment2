package workshop;


import jv.geom.PgBndPolygon;
import jv.geom.PgElementSet;
import jv.geom.PgPolygonSet;
import jv.geom.PgVectorField;
import jv.geom.PuCleanMesh;
import jv.number.PdColor;
import jv.object.PsConfig;
import jv.object.PsDebug;
import jv.object.PsObject;
import jv.project.PgGeometry;
import jv.vecmath.PdVector;
import jv.vecmath.PiVector;
import jv.vecmath.PuMath;
import jv.viewer.PvDisplay;
import jv.project.PvGeometryIf;
import jv.vecmath.PdMatrix;
import jv.vecmath.PdVector;
import jvx.project.PjWorkshop;

import Jama.Matrix;
import Jama.SingularValueDecomposition;

import jvx.project.PjWorkshop;

import java.awt.Color;
import java.util.Vector;

import java.util.Arrays;
import java.util.Collection;
import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.ThreadLocalRandom;

/**
 *  Workshop for surface Interpolation
 */

public class Interpolation extends PjWorkshop {
	
	/** First surface to be registered. */	
	PgElementSet	m_surfP;	
	/** Second surface to be registered. */
	PgElementSet	m_surfQ;	
	
	// Global Variables
	PdVector edge1 = new PdVector(3);
	PdVector edge2 = new PdVector(3);
	PdVector edge3 = new PdVector(3);
	PdVector gx,gy,gz;


	/** Constructor */
	public Interpolation() {
		super("Surface Interpolation");
		if (getClass() == Interpolation.class) {
			init();
		}
	}
	
	/** Initialization */
	public void init() {
		super.init();
	}
	
	// Calculate the normal of triangle surface
	public PdVector calculateTriangleNormal(PdVector v1, PdVector v2, PdVector v3, int ind1, int ind2, int ind3){
		PdVector ab = new PdVector(3);
		PdVector bc = new PdVector(3);
		PdVector normal = new PdVector(3);
		ab = PdVector.subNew(v3,v2);
		edge1 = ab;
		bc = PdVector.subNew(v2,v1);
		edge2 = bc;
		edge3 = PdVector.subNew(v1,v3);
		normal = PdVector.crossNew(ab,bc);
		return normal;
	}
	
	/** Set two Geometries. */
	public void setGeometries(PgElementSet surfP, PgElementSet surfQ) {
		m_surfP = surfP;
		m_surfQ = surfQ;
	}
	
	public void computeMatrixH(int indx1,int indx2,int indx3){
		// Initialize matrices and vectors needed to calculate H
		double A[][] = new double[3][3];
		double B[][] = new double[3][3];
		PdVector v1 = new PdVector(3);
		PdVector v2 = new PdVector(3);
		PdVector v3 = new PdVector(3);
		PdVector w1 = new PdVector(3);
		PdVector w2 = new PdVector(3);
		PdVector w3 = new PdVector(3);
		PdVector vv21 = new PdVector(3);
		PdVector vv31 = new PdVector(3);
		PdVector ww21 = new PdVector(3);
		PdVector ww31 = new PdVector(3);
		
		//Get vertices of the same triangle from input 1 and 2
		v1 = m_surfQ.getVertex(indx1);
		v2 = m_surfQ.getVertex(indx2);
		v3 = m_surfQ.getVertex(indx3);
		w1 = m_surfP.getVertex(indx1);
		w2 = m_surfP.getVertex(indx2);
		w3 = m_surfP.getVertex(indx3);
		vv21 = PdVector.subNew(v2,v1);
		vv31 = PdVector.subNew(v3,v1);
		ww31 = PdVector.subNew(w2,w1);
		ww31 = PdVector.subNew(w3,w1);
		PdVector normal1 = calculateTriangleNormal(v1,v2,v3,indx1,indx2,indx3);
		PdVector normal2 = calculateTriangleNormal(w1,w2,w3,indx1,indx2,indx3);
		A[0][0] =vv21.getEntry(0) ; B[0][0] = ww21.getEntry(0);
		A[0][1] = vv31.getEntry(0); B[0][1] = ww31.getEntry(0);
		A[0][2] = normal1.getEntry(0); B[0][2] = normal2.getEntry(0);
		A[1][0] = vv21.getEntry(1); B[1][0] = ww21.getEntry(1);
		A[1][1] = vv31.getEntry(1); B[1][1] = ww31.getEntry(1);
		A[1][2] = normal1.getEntry(1); B[1][2] = normal2.getEntry(1);
		A[2][0] = vv21.getEntry(2); B[2][0] = ww21.getEntry(2);
		A[2][1] = vv31.getEntry(2); B[2][1] = ww31.getEntry(2);
		A[2][2] = normal1.getEntry(2); B[2][2] = normal2.getEntry(2);

		Matrix matrixA = new Matrix(A,3,3);
		Matrix matrixB = new Matrix(B,3,3);
	}
	
	public void calculationPerTriangle(){
		int numOfElements = m_surfQ.getNumElements();
		int i;
		PiVector indices = new PiVector(3);
		for (i=0;i<numOfElements;i++){
			indices = m_surfQ.getElement(i);
			computeMatrixH(indices.getEntry(0),indices.getEntry(1),indices.getEntry(2));
		}
	}

}

