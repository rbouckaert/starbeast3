package starbeast3.app.beauti;




import java.util.List;


import beast.base.core.BEASTInterface;
import beast.base.core.Function;
import beast.base.core.Input;
import beast.base.inference.StateNode;
import beastfx.app.inputeditor.BeautiDoc;
import beastfx.app.inputeditor.ListInputEditor;
import beastfx.app.util.FXUtils;
import javafx.beans.property.BooleanProperty;
import javafx.beans.property.SimpleBooleanProperty;
import javafx.beans.property.SimpleStringProperty;
import javafx.beans.property.StringProperty;
import javafx.scene.control.CheckBox;
import javafx.scene.control.Label;
import javafx.scene.control.TableColumn;
import javafx.scene.control.TableView;
import javafx.scene.control.cell.CheckBoxTableCell;
import javafx.scene.layout.HBox;
import starbeast3.evolution.branchratemodel.StarBeast3Clock;

public class StarBeast3ClockInputEditor extends ListInputEditor {

	public StarBeast3ClockInputEditor() {
	}
	
	public StarBeast3ClockInputEditor(BeautiDoc doc) {
		super(doc);
	}
	
	@Override
	public Class<?> type() {
		return List.class;
	}

	@Override
	public Class<?> baseType() {
		return StarBeast3Clock.class;
	}
	

	@Override
	public void init(Input<?> input, BEASTInterface beastObject, int itemNr, ExpandOption isExpandOption,
			boolean addButtons) {
		List list = (List)input.get();
		
		pane = FXUtils.newVBox();
		HBox label = FXUtils.newHBox();
		label.getChildren().add(new Label("Select which partition clock rates to estimate. Unestimated rates default to 1. All rates below are relative to the species tree clock rate."));
		pane.getChildren().add(label);
		TableView<Item> table = new TableView<>();


	    CheckBox selecteAllCheckBox = new CheckBox();
	    selecteAllCheckBox.setOnAction(
	        event -> {
	          event.consume();
	          if (doc.autoSetClockRate) {
	        	  return;
	          }
	          table.getItems().forEach(
	        		  item -> item.setSelected(selecteAllCheckBox.isSelected())
	          );
	        });

	    table.setEditable(!doc.autoSetClockRate);
	    
	    TableColumn<Item, String> nameCol = new TableColumn<>("Partition");
	    nameCol.setPrefWidth(550);
	    nameCol.setCellValueFactory(data -> data.getValue().nameProperty());
	    table.getColumns().add(nameCol);
	    
	    TableColumn<Item, Boolean> selectedCol = new TableColumn<>();
	    selectedCol.setGraphic(selecteAllCheckBox);
	    selectedCol.setSortable(false);
	    selectedCol.setPrefWidth(50);

	    selectedCol.setCellFactory( tc -> {
	    	CheckBoxTableCell cell =  new CheckBoxTableCell<>();
	    	cell.setEditable(!doc.autoSetClockRate);
	    	List o = cell.getChildrenUnmodifiable();
	    	return cell;
	    });
	    selectedCol.setCellValueFactory( f -> f.getValue().selectedProperty());
	    selectedCol.setOnEditCommit(e->{
			Boolean newValue = e.getNewValue();
			Item item = (Item) e.getSource();
			item.setSelected(newValue);
		});
	    selectedCol.setEditable(!doc.autoSetClockRate);
	    table.getColumns().add(selectedCol);

	    pane.getChildren().add(table);
	    
		for(Object o : list) {
			StarBeast3Clock clockmodel = (StarBeast3Clock) o;
			table.getItems().add(new Item(clockmodel));
		}
		
	    getChildren().add(pane);
	}

	
	public class Item {
		private StarBeast3Clock clock;
		
	    private BooleanProperty selected = new SimpleBooleanProperty();
	    public void setSelected(Boolean selected) { 
	    	this.selected.set(selected); 
			Function f = clock.meanRateInput.get();
			StateNode s = (StateNode) f;
			s.isEstimatedInput.setValue(selected, s);
	    }
	    public Boolean isSelected() { return selected.get(); }
	    public BooleanProperty selectedProperty() { return selected; }

	    private StringProperty name = new SimpleStringProperty();
	    public void setName(String name) { this.name.set(name); }
	    public String getName() { return name.get(); }
	    public StringProperty nameProperty() { return name; }

	    public Item(StarBeast3Clock clock) {
	    	this.clock = clock;
	    	String id = BeautiDoc.parsePartition(clock.getID());
	    	setName(id);
			Function f = clock.meanRateInput.get();
			if (f instanceof StateNode) {
				StateNode s = (StateNode) f;
				setSelected(s.isEstimatedInput.get());
			}
			selected.addListener((obs,  wasSelected,  isSelected) -> {
				setSelected(isSelected);
			});	
	    }
	  }
}
