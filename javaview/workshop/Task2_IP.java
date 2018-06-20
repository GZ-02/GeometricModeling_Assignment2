package workshop;

import java.awt.*;

import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

import jv.object.PsDialog;
import jv.object.PsUpdateIf;
import jvx.project.PjWorkshop_IP;
import jvx.numeric.PnSparseMatrix;
import jv.vecmath.PdVector;
import jv.vecmath.PdMatrix;

public class Task2_IP extends PjWorkshop_IP implements ActionListener {
	
	//Buttons for assignment 2
	protected Button gradientbtn;
	protected Button laplacebtn;
	protected Button matrixMbtn;
	protected Button matrixSbtn;
	protected Button computeDeformation;
	protected Button resetbtn;
	protected Button testbtn1,testbtn2,testbtn3;

	//Labels for assignment 2
	protected Label gradientLbl;
	protected Label laplaceLbl;
	protected Label matrixMLbl;
	protected Label matrixSLbl;
	protected Label blankLbl,blankLbl1,blankLbl2,blankLbl3,blankLbl4,blankLbl5,blankLbl6,blankLbl7,blankLbl8,blankLbl9,blankLbl10,blankLbl11,blankLbl12;
	protected Label explanationLbl;
	protected Label msgLbl;
	protected Label testLbl1,testLbl2;
	
	// Text fields for assignment 2
	protected TextField txt1,txt2,txt3,txt4,txt5,txt6,txt7,txt8,txt9;

	Task2 t2;
	
	public Task2_IP(){
		super();
		if(getClass()== Task2_IP.class){
			init();
		}
	}
	
	public void init(){
		super.init();
		setTitle("Assignment2");
	}


	public String getNotice() {
		return "Use the mark elements functionality to select triangles of the mesh and then compute deformed surface";
	}
	
	public void setParent(PsUpdateIf parent){
		super.setParent(parent);
		t2 = (Task2)parent;
		
		// Initialize buttons of panel
		gradientbtn = new Button("Calculate Gradient Matrix");
		gradientbtn.addActionListener(this);
		laplacebtn = new Button("Calculate Laplace Matrix");
		laplacebtn.addActionListener(this);
		matrixMbtn = new Button("Calculate Matrix M");
		matrixMbtn.addActionListener(this);
		matrixSbtn = new Button("Calculate Matrix S");
		matrixSbtn.addActionListener(this);
		computeDeformation = new Button ("Compute deformed surface");
		computeDeformation.addActionListener(this);
		resetbtn = new Button("Reset");
		resetbtn.addActionListener(this);
		testbtn1 = new Button("Check G structure");
		testbtn1.addActionListener(this);
		testbtn2 = new Button("Check gradient orthogonal");
		testbtn2.addActionListener(this);

		// Initialize labels of panel
		gradientLbl = new Label();
		laplaceLbl = new Label();
		matrixMLbl = new Label();
		matrixSLbl = new Label();
		explanationLbl = new Label();
		msgLbl = new Label();
		explanationLbl.setText("Specify matrix A:");
		testLbl1 = new Label();
		testLbl2 = new Label();
		// Initialize blank labels that help keeping panel structure
		blankLbl = new Label();
		blankLbl1 =new Label();
		blankLbl2 = new Label();
		blankLbl3 = new Label();
		blankLbl4 = new Label();
		blankLbl5 = new Label();
		blankLbl6 = new Label();
		blankLbl7 = new Label();
		blankLbl8 = new Label();
		blankLbl9 = new Label();
		blankLbl10 = new Label();
		

		// Initialize text fields of panel
		txt1 = new TextField("1.0");
		txt2 = new TextField("0.0");
		txt3 = new TextField("0.0");
		txt4 = new TextField("0.0");
		txt5 = new TextField("1.0");
		txt6 = new TextField("0.0");
		txt7 = new TextField("0.0");
		txt8 = new TextField("0.0");
		txt9 = new TextField("1.0");


		Panel ourPanel = new Panel(new GridLayout(10,4));
		ourPanel.add(gradientbtn);
		ourPanel.add(gradientLbl);
		ourPanel.add(blankLbl);
		ourPanel.add(blankLbl1);
		
		ourPanel.add(laplacebtn);
		ourPanel.add(laplaceLbl);
		ourPanel.add(blankLbl2);
		ourPanel.add(blankLbl3);		
		
		ourPanel.add(matrixMbtn);
		ourPanel.add(matrixMLbl);
		ourPanel.add(blankLbl4);
		ourPanel.add(blankLbl5);

		ourPanel.add(matrixSbtn);
		ourPanel.add(matrixSLbl);
		ourPanel.add(blankLbl6);
		ourPanel.add(blankLbl7);

		ourPanel.add(explanationLbl);
		ourPanel.add(txt1);
		ourPanel.add(txt2);
		ourPanel.add(txt3);

		ourPanel.add(blankLbl8);
		ourPanel.add(txt4);
		ourPanel.add(txt5);
		ourPanel.add(txt6);

		ourPanel.add(blankLbl9);
		ourPanel.add(txt7);
		ourPanel.add(txt8);
		ourPanel.add(txt9);

		ourPanel.add(testbtn1);
		ourPanel.add(testLbl1);
		ourPanel.add(testbtn2);
		ourPanel.add(testLbl2);

		ourPanel.add(computeDeformation);
		ourPanel.add(msgLbl);
		ourPanel.add(blankLbl10);
		ourPanel.add(resetbtn);

		add(ourPanel);

		validate();
		
	}
	
	public void actionPerformed(ActionEvent event) {
		Object source = event.getSource();
		if (source == gradientbtn){
			gradientLbl.setText("...");
			t2.calculateGradientMatrix();
			gradientLbl.setText("Gradient Matrix created");
			return;
		}
		else if (source == laplacebtn){
			laplaceLbl.setText("...");
			t2.calculateinitLaplaceMatrix();
			laplaceLbl.setText("Combinatorial Laplace created");
			return;
		}
		else if (source == matrixMbtn){
			return;
		}
		else if (source == matrixSbtn){
			return;
		}
		else if (source == testbtn1){
			testLbl1.setText("...");
			boolean check1 = t2.TestGradientStructure();
			testLbl1.setText(String.valueOf(check1));
			return;
		}
		else if (source == testbtn2){
			testLbl2.setText("...");
			boolean check2 = t2.TestGradientOrthogonality();
			testLbl2.setText(String.valueOf(check2));
			return;
		}
		else if (source == computeDeformation){
			double a0,a1,a2,a3,a4,a5,a6,a7,a8;
			PdVector result = new PdVector();
			PdVector result2 = new PdVector();
			//double result2;
			a0 = Double.parseDouble(txt1.getText());
			a1 = Double.parseDouble(txt2.getText());
			a2 = Double.parseDouble(txt3.getText());
			a3 = Double.parseDouble(txt4.getText());
			a4 = Double.parseDouble(txt5.getText());
			a5 = Double.parseDouble(txt6.getText());
			a6 = Double.parseDouble(txt7.getText());
			a7 = Double.parseDouble(txt8.getText());
			a8 = Double.parseDouble(txt9.getText());
			if ((a0 == 1.0) && (a1 == 0.0) && (a2 == 0.0) && (a3 == 0.0) && (a4 == 1.0) && (a5 == 0.0) && (a6 == 0.0) && (a7 == 0.0) && (a8 == 1.0)){
				msgLbl.setText("Specify new values for matrix A");
				return;
			}
			double [][] MatrixA = {{a0,a1,a2},{a3,a4,a5},{a6,a7,a8}};
			PdMatrix deformMatrix = new PdMatrix(3, 3);
			deformMatrix.setEntry(0, 0, a0);
			deformMatrix.setEntry(0, 1, a1);
			deformMatrix.setEntry(0, 2, a2);
			deformMatrix.setEntry(1, 0, a3);
			deformMatrix.setEntry(1, 1, a4);
			deformMatrix.setEntry(1, 2, a5);
			deformMatrix.setEntry(2, 0, a6);
			deformMatrix.setEntry(2, 1, a7);
			deformMatrix.setEntry(2, 2, a8);
			msgLbl.setText("...");
			result = t2.calculateGradientEmbeddings();
			result2 = t2.calculateModifiedGradients(MatrixA);
			//compute deformation
			t2.deform(deformMatrix);
			//msgLbl.setText(String.valueOf(result.equals(result2)));
			msgLbl.setText("Done!");
			//t2.m_geom.update(t2.m_geom);
			return;
		}
		else if (source == resetbtn){
			t2.reset();
			return;
		}
	}
	
	protected int getDialogButtons(){
		return PsDialog.BUTTON_OK;
	}
}