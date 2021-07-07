
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;


public class Transcripts {
	
	String transcript_id;
	String parent_gene_id;
	String type;
	String strand; //use "+" and "-" or a boolean
	long gff_mrna_start;
	long bed_mrna_start;
	long mrna_end;
	String mrnachr;
	List<Exon> allExons = new ArrayList<Exon>();
	List<Introns> allIntrons = new ArrayList<Introns>();
	List<CDS> allCDS = new ArrayList<CDS>();
	
	public Transcripts(String tid, String pid, String str, long st, long ed, String c) {
		transcript_id = tid;
		parent_gene_id = pid;
		strand = str;
		gff_mrna_start = st;
		bed_mrna_start = st -1;
		mrna_end = ed;
		mrnachr = c;
	}
	
	public void addExons(Exon exs) {
		allExons.add(exs);
	}
	
	public void addIntrons(Introns ins) {
		allIntrons.add(ins);
	}
	
	public void addCDS(CDS cds) {
		allCDS.add(cds);
	}
	
	public List<Exon> getExons() {
		return allExons;
	}
	
	public List<Introns> getIntrons() {
		return allIntrons;
	}
	
	public List<CDS> getCDS() {
		return allCDS;
	}
}