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
import Jama.EigenvalueDecomposition;
import Jama.Matrix;
import Jama.SingularValueDecomposition;
import jv.geom.PgElementSet;
import jv.object.PsDebug;
import jv.vecmath.PdMatrix;
import jv.vecmath.PdVector;
import jv.vecmath.PiVector;
import jvx.numeric.PnBiconjugateGradient;
import jvx.numeric.PnSparseMatrix;
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
	

	// First model 	
	PgElementSet	meshOrigin;

    
	// Second model
	PgElementSet	meshTarget;
	
    double[] angles;
    PdVector[] rotationAxes;
    PdMatrix[] scalings;
    PdVector[] translations;
    
    public static PdMatrix identity;	
	
	// Global Variables
	PdVector edge1 = new PdVector(3);
	PdVector edge2 = new PdVector(3);
	PdVector edge3 = new PdVector(3);
	PdVector gx,gy,gz;
	PdVector new_gx, new_gy, new_gz;


	/** Constructor */
	public Interpolation() {
		super("Shape Interpolation");
		if (getClass() == Interpolation.class) {
			init();
		}
	}
	
	/** Initialization */
	public void init() {
		super.init();
		identity = new PdMatrix(3);
		identity.setEntry(0, 0, 1);
		identity.setEntry(1, 1, 1);
		identity.setEntry(2, 2, 1);
	}
	

	public void setGeometries(PgElementSet surfP, PgElementSet surfQ) {
		meshOrigin = surfP;
		meshTarget = surfQ;
        calcConstantInfo();
	}
	
	// computation of constant data during interpolation
	public void calcConstantInfo() {
		int numElements = meshOrigin.getNumElements();
		
		angles = new double[numElements];
		rotationAxes = new PdVector[numElements];
		scalings = new PdMatrix[numElements];
		translations = new PdVector[numElements];
		
		PdMatrix[] transforms = getTransforms();
		// SVD for all the elements
		for (int i = 0; i < numElements; i++) {
			SingularValueDecomposition svd = new SingularValueDecomposition(new Matrix(transforms[i].getEntries()));
			
	        Matrix V = svd.getV();
	        Matrix Vt = svd.getV().transpose();
	        Matrix U = svd.getU();
	        Matrix UVt = U.times(Vt);
	        
	        double det = UVt.det();
	        
	        Matrix Z = Matrix.identity(3, 3);
	        Z.set(2, 2, det);
	        
	        Matrix R = U.times(Z).times(Vt);
	        Matrix S = V.times(Z).times(svd.getS()).times(Vt);
	        
	        
	        scalings[i] = new PdMatrix(S.getArray());
	        
	        angles[i] = Math.acos((R.trace()-1)/2.0);
	        
	        Matrix temp = R.plus(R.transpose()).times(0.5);
	        EigenvalueDecomposition eig = temp.eig();
	        for (int v = 0; v < eig.getRealEigenvalues().length; v++){
	        	if (Math.abs(eig.getRealEigenvalues()[v] - 1.0) <= 0.00001 && Math.abs(eig.getImagEigenvalues()[v]) <= 0.00001) {
	        		double[] entries = eig.getV().getColumnPackedCopy();
	        		PdVector axis = new PdVector(3);
	        		axis.setEntry(0, entries[v*3 + 0]);
	        		axis.setEntry(1, entries[v*3 + 1]);
	        		axis.setEntry(2, entries[v*3 + 2]);
	        		rotationAxes[i] = axis;
	        		break;
	        	}
	        }
	        
	        // direction of rotation
	        if (!testAngle(angles[i], rotationAxes[i], R)) {
	        	angles[i] = -angles[i];
	        }
	        
	        PdVector model1 = meshOrigin.getVertex(meshOrigin.getElement(i).getEntry(0));
	        PdVector model2 = meshTarget.getVertex(meshTarget.getElement(i).getEntry(0));
	        
	        translations[i] = PdVector.subNew(model2, model1);
		}
	}
	

	private boolean testAngle(double angle, PdVector axis, Matrix target) {
		PdMatrix temp = getRotationMatrix(angle, axis);
		
		for (int row = 0; row < 3; row++) {
			for (int col = 0; col < 3; col++) {
				double a = temp.getEntry(row, col);
				double b = target.get(row, col);
				if (Math.signum(a) != Math.signum(b)) {
					return false;
				}
			}
		}
		
		return true;
	}
	
	//this creates the looseMesh from time argument
	public PgElementSet getInterpolatedset(double time) {
		int numElements = meshOrigin.getNumElements();
		
		// Initiate the new geometry
		PgElementSet newSet = new PgElementSet();
		newSet.setNumElements(numElements);
		newSet.setNumVertices(numElements*3);
		
		PdVector[] normals = new PdVector[numElements];
		
		for (int index = 0; index < numElements; index++) {
			PiVector faceOrigin = meshOrigin.getElement(index);
			PdVector p1 = meshOrigin.getVertex(faceOrigin.getEntry(0));
			PdVector p2 = meshOrigin.getVertex(faceOrigin.getEntry(1));
			PdVector p3 = meshOrigin.getVertex(faceOrigin.getEntry(2));
			PdVector nv = meshOrigin.getElementNormal(index);
			
			PdVector v1 = PdVector.subNew(p2, p1);
			PdVector v2 = PdVector.subNew(p3, p1);
			
			PdMatrix scaler = new PdMatrix(3);
			scaler.multScalar(scalings[index], time);

			PdMatrix id = new PdMatrix(3);
			id.multScalar(identity, 1-time);
			
			id.add(scaler);
			

			PdMatrix rot = getRotationMatrix(time * angles[index], rotationAxes[index]);

			rot.rightMult(id);
			
			v1.leftMultMatrix(rot);
			v2.leftMultMatrix(rot);
			nv.leftMultMatrix(rot);

			v1.add(p1);
			v2.add(p1);
			
			PdVector translation = new PdVector(translations[index].getEntries());
			translation.multScalar(time);

			PdVector p1New = PdVector.addNew(p1, translation);
			PdVector p2New = PdVector.addNew(v1, translation);
			PdVector p3New = PdVector.addNew(v2, translation);
			
			newSet.setVertex(index*3 + 0, p1New);
			newSet.setVertex(index*3 + 1, p2New);
			newSet.setVertex(index*3 + 2, p3New);

			newSet.setElement(index, new PiVector(index*3, index*3+1, index*3+2));
			
			nv.normalize();
			normals[index] = nv;
		}
		
		newSet.setElementNormals(normals);
		
		newSet.update(newSet);
		
		return newSet;
	}
	

	private PdMatrix getRotationMatrix(double angle, PdVector axis) {
		double C = Math.cos(angle);
		double S = Math.sin(angle);
		
		double Cinv = 1.0 - C;
		
		double a1 = axis.getEntry(0);
		double a2 = axis.getEntry(1);
		double a3 = axis.getEntry(2);
		
		PdMatrix result = new PdMatrix(3);
		result.setEntry(0, 0, Math.pow(a1, 2) + C*(1-Math.pow(a1, 2)));
		result.setEntry(1, 0, a1*a2*Cinv + a3*S);
		result.setEntry(2, 0, a1*a3*Cinv - a2*S);
		
		result.setEntry(0, 1, a1*a2*Cinv - a3*S);
		result.setEntry(1, 1, Math.pow(a2, 2)+ C*(1-Math.pow(a2, 2)));
		result.setEntry(2, 1, a2*a3*Cinv + a1*S);
		
		result.setEntry(0, 2, a1*a3*Cinv + a2*S);
		result.setEntry(1, 2, a2*a3*Cinv - a1*S);
		result.setEntry(2, 2, Math.pow(a3, 2) + C*(1-Math.pow(a3, 2)));
		
		return result;
	}
	// H matrix
	private PdMatrix[] getTransforms() {
		meshOrigin.makeElementNormals();
		meshTarget.makeElementNormals();
		
		PdMatrix[] transforms = new PdMatrix[meshOrigin.getNumElements()];
		
		for (int i = 0; i < meshOrigin.getNumElements(); i++) {
			PdMatrix A = getTransform(meshOrigin, i);
			PdMatrix B = getTransform(meshTarget, i);
			
			A.invert();
			B.rightMult(A);
			
			
			transforms[i] = B;
		}
		
		return transforms;
	}
	
	private PdMatrix getTransform(PgElementSet mesh, int index) {
		PiVector faceOrigin = mesh.getElement(index);
		PdVector p1 = mesh.getVertex(faceOrigin.getEntry(0));
		PdVector p2 = mesh.getVertex(faceOrigin.getEntry(1));
		PdVector p3 = mesh.getVertex(faceOrigin.getEntry(2));
		PdVector nv = mesh.getElementNormal(index);
		
		PdVector v1 = PdVector.subNew(p2, p1);
		PdVector v2 = PdVector.subNew(p3, p1);
		
		PdMatrix V = new PdMatrix(3, 3);
		V.setColumn(0, v1);
		V.setColumn(1, v2);
		V.setColumn(2, nv);
		
		return V;
	}
	//called from interface
	public PgElementSet getGradientInterpolated(PgElementSet looseMesh) {
		return interpolateSet(meshOrigin, looseMesh);
	}
	
	private PgElementSet interpolateSet(PgElementSet origin, PgElementSet intermediate) {
		PgElementSet copy = (PgElementSet) origin.clone();

        PnSparseMatrix G = meshToGradient(origin);
		PnSparseMatrix GT = PnSparseMatrix.transposeNew(G);
        PnSparseMatrix M = calculateM();
    	PnSparseMatrix left = PnSparseMatrix.multMatrices(GT, PnSparseMatrix.multMatrices(M, G, null), null);
    	PnSparseMatrix leftHand = PnSparseMatrix.copyNew(left);

    	PdVector x = new PdVector(origin.getNumVertices());
    	PdVector y = new PdVector(origin.getNumVertices());
    	PdVector z = new PdVector(origin.getNumVertices());
        
        PdVector[] g = meshToGradientVector(origin, intermediate);
        PnSparseMatrix right = PnSparseMatrix.multMatrices(GT, M, null);
        
        PdVector xGradient = PnSparseMatrix.rightMultVector(right, g[0], null);
        PdVector yGradient = PnSparseMatrix.rightMultVector(right, g[1], null);
        PdVector zGradient = PnSparseMatrix.rightMultVector(right, g[2], null);

    	PnBiconjugateGradient CGSOLVER = new PnBiconjugateGradient();
    		
    	CGSOLVER.solve(leftHand, x, xGradient);
    	CGSOLVER.solve(leftHand, y, yGradient);
    	CGSOLVER.solve(leftHand, z, zGradient);

    	PdVector sumNew = new PdVector(3);
    	for (int vIndex = 0; vIndex < origin.getNumVertices(); vIndex++) {
    		sumNew.setEntry(0, sumNew.getEntry(0) + x.getEntry(vIndex));
    		sumNew.setEntry(1, sumNew.getEntry(1) + y.getEntry(vIndex));
    		sumNew.setEntry(2, sumNew.getEntry(2) + z.getEntry(vIndex));
    	}
    	sumNew.multScalar(1.0 / origin.getNumVertices());
    	PdVector sumOld = getMean(intermediate.getVertices());
    	
    	PdVector translationMean = PdVector.subNew(sumOld, sumNew);
    	
    	for (int vIndex = 0; vIndex < copy.getNumVertices(); vIndex++) {
    		PdVector newV = new PdVector(3);
    		newV.setEntry(0, x.getEntry(vIndex));
    		newV.setEntry(1, y.getEntry(vIndex));
    		newV.setEntry(2, z.getEntry(vIndex));
    		
    		newV.add(translationMean);
    		
    		copy.setVertex(vIndex, newV);
    	}
    	
		return copy;
	}
	
	private PdVector getMean(PdVector[] vertices) {
		PdVector sum = new PdVector(3);
		for (PdVector v : vertices) {
			sum.add(v);
		}
		sum.multScalar(1.0 / vertices.length);
		return sum;
	}
	
	
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
		PdVector cross1 = new PdVector();
		PdVector cross2 = new PdVector();
		PdVector cross3 = new PdVector();

		// Get vertices of a triangle
		v1 = meshOrigin.getVertex(index1);
		v2 = meshOrigin.getVertex(index2);
		v3 = meshOrigin.getVertex(index3);

		PdVector normal = calculateTriangleNormal(v1,v2,v3,index1,index2,index3);
		double embadon = calculateTriangleArea(v1,v2,v3);

		cross1 = PdVector.crossNew(normal,edge1);
		cross2 = PdVector.crossNew(normal,edge2);
		cross3 = PdVector.crossNew(normal,edge3); 

		triangleGradient[0][0] = (1/(2*embadon)) * cross1.getEntry(0); 
		triangleGradient[0][1] = (1/(2*embadon)) * cross2.getEntry(0); 
		triangleGradient[0][2] = (1/(2*embadon)) * cross3.getEntry(0); 
		triangleGradient[1][0] = (1/(2*embadon)) * cross1.getEntry(1); 
		triangleGradient[1][1] = (1/(2*embadon)) * cross2.getEntry(1); 
		triangleGradient[1][2] = (1/(2*embadon)) * cross3.getEntry(1); 
		triangleGradient[2][0] = (1/(2*embadon)) * cross1.getEntry(2); 
		triangleGradient[2][1] = (1/(2*embadon)) * cross2.getEntry(2); 
		triangleGradient[2][2] = (1/(2*embadon)) * cross3.getEntry(2); 
		return triangleGradient;
	}

	public PnSparseMatrix calculateGradientMatrix(){
		// Declaration of variables
		int numOfElements = meshOrigin.getNumElements();
		int numOfVertices = meshOrigin.getNumVertices();
		PiVector VertexIndices = new PiVector(3);
		int index1,index2,index3;
		double gradient[][] = new double[3][3];
		PnSparseMatrix gradientMatrix = new PnSparseMatrix(3*numOfElements,numOfVertices,3);

		for (int i = 0; i<numOfElements;i++){
			// Find indices of triangle vertices
			VertexIndices = meshOrigin.getElement(i);
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
	
	public PnSparseMatrix calculateM() {
    	
		PnSparseMatrix M = new PnSparseMatrix(meshOrigin.getNumElements() * 3, meshOrigin.getNumElements() * 3, 1);
    	PiVector [] triangleArray = meshOrigin.getElements();
    	
		for(int triangleIndex = 0; triangleIndex < triangleArray.length; triangleIndex++) {
            PiVector triangle = triangleArray[triangleIndex];
            
            PdVector v1 = meshOrigin.getVertex(triangle.getEntry(0));
            PdVector v2 = meshOrigin.getVertex(triangle.getEntry(1));
            PdVector v3 = meshOrigin.getVertex(triangle.getEntry(2));

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
		PnSparseMatrix M = calculateM();
        PnSparseMatrix GT = PnSparseMatrix.transposeNew(G);
        
        // S=G^T*M*G
		PnSparseMatrix S = PnSparseMatrix.multMatrices(M, G, null);
		S = PnSparseMatrix.multMatrices(GT, S, null);
		
		return S;
		
	}
	
		public static PdVector[] meshToGradientVector(PgElementSet mesh, PgElementSet gradientTarget) {
        
        PiVector[] trianglesMesh = mesh.getElements();
        PiVector[] trianglesTarget = gradientTarget.getElements();
        
        PdVector x = new PdVector(mesh.getNumElements()*3);
        PdVector y = new PdVector(mesh.getNumElements()*3);
        PdVector z = new PdVector(mesh.getNumElements()*3);

        for(int triangleIdx = 0; triangleIdx < trianglesMesh.length; triangleIdx++) {
            PiVector triangleMesh = trianglesMesh[triangleIdx];
            PiVector triangleTarget = trianglesTarget[triangleIdx];
            
            PdVector[] vertices = new PdVector[]{
            		mesh.getVertex(triangleMesh.getEntry(0)),
            		mesh.getVertex(triangleMesh.getEntry(1)),
            		mesh.getVertex(triangleMesh.getEntry(2))};

            PdMatrix subGradient = triangleToGradient(vertices,
            		mesh.getElementNormal(triangleIdx));
            
            PdVector xTemp = new PdVector(3);
            PdVector yTemp = new PdVector(3);
            PdVector zTemp = new PdVector(3);
            
            for (int i = 0; i < vertices.length; i++) {
            	PdVector v = gradientTarget.getVertex(triangleTarget.getEntry(i));
            	xTemp.setEntry(i, v.getEntry(0));
            	yTemp.setEntry(i, v.getEntry(1));
            	zTemp.setEntry(i, v.getEntry(2));
            }
            
            xTemp.leftMultMatrix(subGradient);
            yTemp.leftMultMatrix(subGradient);
            zTemp.leftMultMatrix(subGradient);
            
            for (int i = 0; i < 3; i++) {
            	x.setEntry(triangleIdx * 3 + i, xTemp.getEntry(i));
            	y.setEntry(triangleIdx * 3 + i, yTemp.getEntry(i));
            	z.setEntry(triangleIdx * 3 + i, zTemp.getEntry(i));
            }
        }

        return new PdVector[] {x, y, z};
    }
	
	private static PdMatrix triangleToGradient(PdVector[] vertices, PdVector normal) {
        PdMatrix gradient = new PdMatrix(3, 3);

        PdVector p1 = vertices[0];
        PdVector p2 = vertices[1];
        PdVector p3 = vertices[2];

        PdVector V = PdVector.subNew(p2, p1);
        PdVector W = PdVector.subNew(p3, p1);

       
        double area = PdVector.crossNew(V, W).length() * 0.5;

        PdVector e1 = PdVector.subNew(p3, p2);
        PdVector e2 = PdVector.subNew(p1, p3);
        PdVector e3 = PdVector.subNew(p2, p1);

        gradient.setColumn(0, PdVector.crossNew(normal, e1));
        gradient.setColumn(1, PdVector.crossNew(normal, e2));
        gradient.setColumn(2, PdVector.crossNew(normal, e3));

        gradient.multScalar(1.0 / (area * 2));

        return gradient;
    }
	
	public static PnSparseMatrix meshToGradient(PgElementSet mesh) {
        PnSparseMatrix G = new PnSparseMatrix(mesh.getNumElements() * 3, mesh.getNumVertices(), 3);
        PiVector[] triangles = mesh.getElements();

        for(int triangleIdx = 0; triangleIdx < triangles.length; triangleIdx++) {
            PiVector triangle = triangles[triangleIdx];

            PdMatrix subGradient = triangleToGradient(new PdVector[]{
            		mesh.getVertex(triangle.getEntry(0)),
            		mesh.getVertex(triangle.getEntry(1)),
            		mesh.getVertex(triangle.getEntry(2))},
            		mesh.getElementNormal(triangleIdx));

            for(int columnIdx = 0; columnIdx < 3; columnIdx++) {
                int column = 3 * triangleIdx;

                G.addEntry(column, triangle.getEntry(columnIdx), subGradient.getColumn(columnIdx).getEntry(0));
                G.addEntry(column + 1, triangle.getEntry(columnIdx), subGradient.getColumn(columnIdx).getEntry(1));
                G.addEntry(column + 2, triangle.getEntry(columnIdx), subGradient.getColumn(columnIdx).getEntry(2));
            }
        }

        return G;
    }

}

