package potts.likelihood;

import java.io.PrintStream;

import org.w3c.dom.Node;

import beast.base.inference.StateNode;

public class SequenceState extends StateNode {

	private int [][][] sequences;
	private int [] currentSequence;
	private int [] storedSequence;
	

	public void init(int nodeCount, int siteCount) {
		sequences = new int[2][nodeCount][siteCount];
		currentSequence = new int [nodeCount];
		storedSequence = new int [nodeCount];
	}
	
	public void setSequence(int i, int [] states) {
        System.arraycopy(states, 0, sequences[0][i], 0, states.length);
	}
	
	public int [] get(int i) {
		return sequences[currentSequence[i]][i];
	}

	@Override
	public void init(PrintStream out) {
		// TODO Auto-generated method stub

	}

	@Override
	public void log(long sample, PrintStream out) {
		// TODO Auto-generated method stub

	}

	@Override
	public void close(PrintStream out) {
		// TODO Auto-generated method stub

	}

	@Override
	public int getDimension() {
		// TODO Auto-generated method stub
		return 0;
	}

	@Override
	public double getArrayValue(int dim) {
		// TODO Auto-generated method stub
		return 0;
	}

	@Override
	public void initAndValidate() {
		// TODO Auto-generated method stub

	}

	@Override
	public void setEverythingDirty(boolean isDirty) {
		// TODO Auto-generated method stub

	}

	@Override
	public StateNode copy() {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public void assignTo(StateNode other) {
		// TODO Auto-generated method stub

	}

	@Override
	public void assignFrom(StateNode other) {
		// TODO Auto-generated method stub

	}

	@Override
	public void assignFromFragile(StateNode other) {
		// TODO Auto-generated method stub

	}

	@Override
	public void fromXML(Node node) {
		// TODO Auto-generated method stub

	}

	@Override
	public int scale(double scale) {
		// TODO Auto-generated method stub
		return 0;
	}

	@Override
	protected void store() {
		System.arraycopy(currentSequence, 0, storedSequence, 0, currentSequence.length);
	}

	@Override
	public void restore() {
		int [] tmp = currentSequence;
		currentSequence = storedSequence;
		storedSequence = tmp;
	}

}
