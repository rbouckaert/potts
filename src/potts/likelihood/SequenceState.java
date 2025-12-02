package potts.likelihood;

import java.io.PrintStream;

import org.w3c.dom.Node;

import beast.base.core.Description;
import beast.base.inference.StateNode;

@Description("State node representing sequences on all nodes (including internal nodes) on a tree")
public class SequenceState extends StateNode {

	private int [][][] sequences;
	private int [] currentSequenceIndex;
	private int [] storedSequenceIndex;
	

	public void init(int nodeCount, int siteCount) {
		sequences = new int[2][nodeCount][siteCount];
		currentSequenceIndex = new int [nodeCount];
		storedSequenceIndex = new int [nodeCount];
	}
	
	public boolean isDirtySequence(int i) {
		return storedSequenceIndex[i] != currentSequenceIndex[i];
	}
	
	public void setSequence(int i, int [] states) {
        System.arraycopy(states, 0, sequences[0][i], 0, states.length);
	}
	
	public void setEditedSequence(int i, int [] states) {
		startEditing(null);
		currentSequenceIndex[i] = 1-currentSequenceIndex[i];
        System.arraycopy(states, 0, sequences[currentSequenceIndex[i]][i], 0, states.length);
	}

	
	public int [] get(int i) {
		return sequences[currentSequenceIndex[i]][i];
	}


	@Override
	protected void store() {
		System.arraycopy(currentSequenceIndex, 0, storedSequenceIndex, 0, currentSequenceIndex.length);
	}

	@Override
	public void restore() {
		int [] tmp = currentSequenceIndex;
		currentSequenceIndex = storedSequenceIndex;
		storedSequenceIndex = tmp;
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
}
