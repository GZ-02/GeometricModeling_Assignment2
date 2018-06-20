package workshop;

import java.util.Stack;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Set;
import jv.geom.PgEdgeStar;
import jv.geom.PgElementSet;
import jv.object.PsDebug;
import jv.project.PgGeometry;
import jv.vecmath.PdVector;
import jv.vecmath.PiVector;
import jvx.project.PjWorkshop;

import jvx.numeric.PnSparseMatrix;
import java.lang.Math;
import jv.object.PsObject;
import jv.geom.PgPolygon;

public class Task2 extends PjWorkshop {
	PgElementSet m_geom;
	PgElementSet m_geomSave;
	boolean traversedAll = false; //for checking if we traversed all the vertices in the geometry
	
	// Global Variables
	PdVector edge1 = new PdVector(3);
	PdVector edge2 = new PdVector(3);
	PdVector edge3 = new PdVector(3);
	PdVector gx,gy,gz;
	PdVector new_gx, new_gy, new_gz; 

	public Task2() {
		super("Task2 of the geometric model");
		init();
	}
	
	@Override
	public void setGeometry(PgGeometry geom) {
		super.setGeometry(geom);
		m_geom 		= (PgElementSet)super.m_geom;
		m_geomSave 	= (PgElementSet)super.m_geomSave;
	}
	
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

	public double calculateTriangleArea(PdVector v1, PdVector v2, PdVector v3){
		// Calculate area of triangle using Heron's formula
		// https://www.mathopenref.com/heronsformula.html
		double a,b,c;
		a = Math.pow(v1.getEntry(0)-v2.getEntry(0),2) + Math.pow(v1.getEntry(1)-v2.getEntry(1),2) + Math.pow(v1.getEntry(2)-v2.getEntry(2),2);
		a = Math.sqrt(a);
		b = Math.pow(v1.getEntry(0)-v3.getEntry(0),2) + Math.pow(v1.getEntry(1)-v3.getEntry(1),2) + Math.pow(v1.getEntry(2)-v3.getEntry(2),2);
		b = Math.sqrt(b);
		c = Math.pow(v2.getEntry(0)-v3.getEntry(0),2) + Math.pow(v2.getEntry(1)-v3.getEntry(1),2) + Math.pow(v2.getEntry(2)-v3.getEntry(2),2);
		c = Math.sqrt(c);
		double p = (a+b+c)/2;
		double area  = p*(p-a)*(p-b)*(p-c); // calculate half of the perimeter of the triangle
		area = Math.sqrt(area);
		return area;
	}

	public double[][] calculateGradientOfTriangle(int index1,int index2,int index3){
		double triangleGradient [][] = new double [3][3];
		// Initialize vectors showing vertices of triangle
		PdVector v1 = new PdVector(3);
		PdVector v2 = new PdVector(3);
		PdVector v3 = new PdVector(3);
		
		// Get vertices of a triangle
		v1 = m_geom.getVertex(index1);
		v2 = m_geom.getVertex(index2);
		v3 = m_geom.getVertex(index3);

		PdVector normal = calculateTriangleNormal(v1,v2,v3,index1,index2,index3);
		double embadon = calculateTriangleArea(v1,v2,v3);
		triangleGradient[0][0] = (1/(2*embadon)) * normal.getEntry(0)*edge1.getEntry(0);
		triangleGradient[0][1] = (1/(2*embadon)) * normal.getEntry(0)*edge2.getEntry(0);
		triangleGradient[0][2] = (1/(2*embadon)) * normal.getEntry(0)*edge3.getEntry(0);
		triangleGradient[1][0] = (1/(2*embadon)) * normal.getEntry(1)*edge1.getEntry(1);
		triangleGradient[1][1] = (1/(2*embadon)) * normal.getEntry(1)*edge2.getEntry(1);
		triangleGradient[1][2] = (1/(2*embadon)) * normal.getEntry(1)*edge3.getEntry(1);
		triangleGradient[2][0] = (1/(2*embadon)) * normal.getEntry(2)*edge1.getEntry(2);
		triangleGradient[2][1] = (1/(2*embadon)) * normal.getEntry(2)*edge2.getEntry(2);
		triangleGradient[2][2] = (1/(2*embadon)) * normal.getEntry(2)*edge3.getEntry(2);
		return triangleGradient;
	}

	public PnSparseMatrix calculateGradientMatrix(){
		// Declaration of variables
		int numOfElements = m_geom.getNumElements();
		int numOfVertices = m_geom.getNumVertices();
		PiVector VertexIndices = new PiVector(3);
		int index1,index2,index3;
		double gradient[][] = new double[3][3];
		PnSparseMatrix gradientMatrix = new PnSparseMatrix(3*numOfElements,numOfVertices,3);

		for (int i = 0; i<numOfElements;i++){
			// Find indices of triangle vertices
			VertexIndices = m_geom.getElement(i);
			index1 = VertexIndices.getEntry(0);
			index2 = VertexIndices.getEntry(1);
			index3 = VertexIndices.getEntry(2);
			gradient = calculateGradientOfTriangle(index1,index2,index3);
			// Create gradient matrix for mesh
			// Three entries in first row of triangle i
			gradientMatrix.addEntry(3*i+0,index1,gradient[0][0]);
			gradientMatrix.addEntry(3*i+0,index2,gradient[0][1]);
			gradientMatrix.addEntry(3*i+0,index3,gradient[0][2]);
			// Three entries in second row of triangle i
			gradientMatrix.addEntry(3*i+1,index1,gradient[1][0]);
			gradientMatrix.addEntry(3*i+1,index2,gradient[1][1]);
			gradientMatrix.addEntry(3*i+1,index3,gradient[1][2]);
			// Three entries in third row of triangle i
			gradientMatrix.addEntry(3*i+2,index1,gradient[2][0]);
			gradientMatrix.addEntry(3*i+2,index2,gradient[2][1]);
			gradientMatrix.addEntry(3*i+2,index3,gradient[2][2]);
		}
		return gradientMatrix;	
	}

	public boolean TestGradientStructure(){
		boolean test = false;
		PnSparseMatrix G = calculateGradientMatrix();
		PiVector numOfnonZeros = new PiVector(); 
		numOfnonZeros = G.getNumEntries();
		int i;
		for(i=0;i<numOfnonZeros.getSize();i++){
			if (numOfnonZeros.getEntry(i) == 3)
				test = true;
			else
				test=false;
		}
		return test;
	}

	public boolean TestGradientOrthogonality(){
		boolean test = false;
		boolean selected;
		int numOfElements = m_geom.getNumElements();
		PiVector VertexIndices = new PiVector(3);
		int index1,index2,index3,i;
		double gradient[][] = new double[3][3];
		int count=0;
		for(i=0;i<numOfElements;i++){
			VertexIndices = m_geom.getElement(i);
			selected = VertexIndices.hasTag(PsObject.IS_SELECTED);
			if(selected){
				VertexIndices = m_geom.getElement(i);
				index1 = VertexIndices.getEntry(0);
				index2 = VertexIndices.getEntry(1);
				index3 = VertexIndices.getEntry(2);
				gradient = calculateGradientOfTriangle(index1,index2,index3);
				PdVector grad1 = new PdVector(gradient[0][0],gradient[1][0],gradient[2][0]);
				PdVector grad2 = new PdVector(gradient[0][1],gradient[1][1],gradient[2][1]);
				PdVector grad3 = new PdVector(gradient[0][2],gradient[1][2],gradient[2][2]);
				if ((grad1.dot(edge1) == 0) && (grad2.dot(edge2) == 0) && (grad3.dot(edge3) == 0)){
					test=true;
					count++;
				}
				else{
					test=false;	
				} 
			}
		}
		return test;
	}


	// Calculates sparse matrix L in the formula D = L.v(x,y,z) {Lecture 5, slide 8}
	// Also refer figure 3 in https://people.eecs.berkeley.edu/~jrs/meshpapers/Sorkine.pdf
	public PnSparseMatrix calculateinitLaplaceMatrix(){
		
		int numOfVertices = m_geom.getNumVertices();
		PnSparseMatrix initLaplace = new PnSparseMatrix(numOfVertices,numOfVertices);
		//PdVector [] vertices = m_geom.getVertices();
		PiVector valence = new PiVector();
		valence = m_geom.getVertexValence(m_geom);
		for(int i = 0; i < numOfVertices; i++)
		{	
			PiVector neighbours = m_geom.getNeighbour(i);
			initLaplace.addEntry(i,i,valence.getEntry(i));
			for(int j =0; j<neighbours.getSize(); j++){
				if(initLaplace.getEntry(i,neighbours.getEntry(j)) == -1)
					continue;
				initLaplace.addEntry(i,neighbours.getEntry(j),-1);
				initLaplace.addEntry(neighbours.getEntry(j),i,-1);
			}	
		}
		return initLaplace;
	}

	// Calculate D=L.v
	public PdVector Laplacemult(PnSparseMatrix L, PdVector v){
		PnSparseMatrix result = new PnSparseMatrix(m_geom.getNumVertices(),m_geom.getNumVertices());
		return result.rightMultVector(L, v, null);
	}

	// Calculate combinatorial Laplace
	// Have to define 3 vertex arrays for x,y,z
	public void CombiLaplace(PdVector a, PdVector b, PdVector c){

		PdVector dx = new PdVector(a.getSize());
		PdVector dy = new PdVector(b.getSize());
		PdVector dz = new PdVector(c.getSize());

		PnSparseMatrix L = calculateinitLaplaceMatrix();

		dx = Laplacemult(L, a);
		dy = Laplacemult(L, b);
		dz = Laplacemult(L, c);
	}
	
	
	public PnSparseMatrix calculateMatrixM() {
    	
		PnSparseMatrix M = new PnSparseMatrix(m_geom.getNumElements() * 3, m_geom.getNumElements() * 3, 1);
    	PiVector [] triangleArray = m_geom.getElements();
    	
		for(int triangleIndex = 0; triangleIndex < triangleArray.length; triangleIndex++) {
            PiVector triangle = triangleArray[triangleIndex];
            
            PdVector v1 = m_geom.getVertex(triangle.getEntry(0));
            PdVector v2 = m_geom.getVertex(triangle.getEntry(1));
            PdVector v3 = m_geom.getVertex(triangle.getEntry(2));

            PdVector A = PdVector.subNew(v2, v1);
            PdVector B = PdVector.subNew(v3, v1);
 
            double area = PdVector.crossNew(A, B).length() * 0.5; // area = ||(p2 - p1) x (p3 - p1)|| * 0.5
            //area = area/3; //Not sure if one third of the area is needed
            int position = triangleIndex * 3;
            
            M.setEntry(position, position, area);
            M.setEntry(position+1, position+1, area);
            M.setEntry(position+2, position+2, area);
    	}
    	return M;
    }
	
	public PnSparseMatrix calculateMatrixS(){
		PnSparseMatrix G = calculateGradientMatrix();
		PnSparseMatrix M = calculateMatrixM();
        PnSparseMatrix GT = PnSparseMatrix.transposeNew(G);
        
        // S=G^T*M*G
		PnSparseMatrix S = PnSparseMatrix.multMatrices(M, G, null);
		S = PnSparseMatrix.multMatrices(GT, S, null);
		
		return S;
		
	}

	//Function that calculates the gradient emebeddings for the mesh
	public PdVector calculateGradientEmbeddings(){
		int numOfElements = m_geom.getNumElements();
		int numOfVertices = m_geom.getNumVertices();
		PnSparseMatrix G = calculateGradientMatrix();
		PdVector vx = new PdVector(numOfVertices);
		PdVector vy = new PdVector(numOfVertices);
		PdVector vz = new PdVector(numOfVertices);
		gx = new PdVector(3*numOfElements);
		gy = new PdVector(3*numOfElements);
		gz = new PdVector(3*numOfElements);
		PdVector vertex = new PdVector(3);
		int i;
		for(i=0;i<numOfVertices;i++){
			vertex = m_geom.getVertex(i);
			vx.setEntry(i,vertex.getEntry(0));
			vy.setEntry(i,vertex.getEntry(1));
			vz.setEntry(i,vertex.getEntry(2));
		}
		PnSparseMatrix.rightMultVector(G,vx,gx);
		PnSparseMatrix.rightMultVector(G,vy,gy);
		PnSparseMatrix.rightMultVector(G,vz,gz);
		return gx;
	}


	// Calculate modified gradients
	public PdVector calculateModifiedGradients(double[][] A){
		boolean selected;
		int numOfElements = m_geom.getNumElements();
		int i;		
		double x1,x2,x3;
		double y1,y2,y3;
		double z1,z2,z3;
		new_gx = new PdVector(3*numOfElements);
		new_gy = new PdVector(3*numOfElements);
		new_gz = new PdVector(3*numOfElements);
		for (i=0;i<3*numOfElements;i++){
			new_gx.setEntry(i,gx.getEntry(i));
			new_gy.setEntry(i,gy.getEntry(i));
			new_gz.setEntry(i,gz.getEntry(i));
		}
		for(i=0;i<numOfElements;i++){
			PiVector element = m_geom.getElement(i);
			selected = element.hasTag(PsObject.IS_SELECTED);
			// Check if the triangle was selected by the user
			// If it was, modify the appropriate gradients 
			if(selected){
				// Modify gradient for x coordinates
				x1 = A[0][0]*gx.getEntry(3*i+0) + A[0][1]*gx.getEntry(3*i+1) + A[0][2]*gx.getEntry(3*i+2);
				x2 = A[1][0]*gx.getEntry(3*i+0) + A[1][1]*gx.getEntry(3*i+1) + A[1][2]*gx.getEntry(3*i+2);
				x3 = A[2][0]*gx.getEntry(3*i+0) + A[2][1]*gx.getEntry(3*i+1) + A[2][2]*gx.getEntry(3*i+2);
				new_gx.setEntry(3*i+0,x1);
				new_gx.setEntry(3*i+1,x2);
				new_gx.setEntry(3*i+2,x3);

				// Modify gradients for y coordinates
				y1 = A[0][0]*gy.getEntry(3*i+0) + A[0][1]*gy.getEntry(3*i+1) + A[0][2]*gy.getEntry(3*i+2);
				y2 = A[1][0]*gy.getEntry(3*i+0) + A[1][1]*gy.getEntry(3*i+1) + A[1][2]*gy.getEntry(3*i+2);
				y3 = A[2][0]*gy.getEntry(3*i+0) + A[2][1]*gy.getEntry(3*i+1) + A[2][2]*gy.getEntry(3*i+2);
				new_gy.setEntry(3*i+0,y1);
				new_gy.setEntry(3*i+1,y2);
				new_gy.setEntry(3*i+2,y3);
				
				// Modify gradients for z coordinates
				z1 = A[0][0]*gz.getEntry(3*i+0) + A[0][1]*gz.getEntry(3*i+1) + A[0][2]*gz.getEntry(3*i+2);
				z2 = A[1][0]*gz.getEntry(3*i+0) + A[1][1]*gz.getEntry(3*i+1) + A[1][2]*gz.getEntry(3*i+2);
				z3 = A[2][0]*gz.getEntry(3*i+0) + A[2][1]*gz.getEntry(3*i+1) + A[2][2]*gz.getEntry(3*i+2);
				new_gz.setEntry(3*i+0,z1);
				new_gz.setEntry(3*i+1,z2);
				new_gz.setEntry(3*i+2,z3);
			}
		}
		return new_gx;
	}


	// method that implements simplified version of shape editing
	public void computeDeformation(double[][] A){
		PdVector rightSideX = new PdVector();
		PdVector rightSideY = new PdVector();
		PdVector rightSideZ = new PdVector();
		PdVector newVx = new PdVector();
		PdVector newVy = new PdVector();
		PdVector newVz = new PdVector();
		// First we calculate the gradients of embedding gx,gy and gz
		calculateGradientEmbeddings();		
		// Then based on the selection of user we modify the gradients and get new ones ~gx,~gy and ~gz
		calculateModifiedGradients(A);		
		// Get matrices M,S and G to 
		PnSparseMatrix S = calculateMatrixS();
		PnSparseMatrix M = calculateMatrixM();
		PnSparseMatrix G = calculateGradientMatrix();
		PnSparseMatrix GT = PnSparseMatrix.transposeNew(G);

		PnSparseMatrix B = PnSparseMatrix.multMatrices(GT,M,null);
		PnSparseMatrix.rightMultVector(B,new_gx,rightSideX);
		PnSparseMatrix.rightMultVector(B,new_gy,rightSideY);
		PnSparseMatrix.rightMultVector(B,new_gz,rightSideZ);
}

}
