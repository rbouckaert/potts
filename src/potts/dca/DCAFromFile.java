package potts.dca;

import java.io.File;
import java.io.IOException;

import org.json.JSONException;
import org.json.JSONObject;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beastfx.app.inputeditor.BeautiDoc;

@Description("Load DCA that was previously trained from file")
public class DCAFromFile extends DCA {

	final public Input<File> dcaInput = new Input<>("fileName", "DCA json file previously trained on an alignment",
			Validate.REQUIRED);

	@Override
	public void initAndValidate() {
		String json;
		try {
			json = BeautiDoc.load(dcaInput.get());
			fromJSON(new JSONObject(json));
		} catch (IOException | JSONException e) {
			e.printStackTrace();
		}
		
		super.initAndValidate();
	}
}
