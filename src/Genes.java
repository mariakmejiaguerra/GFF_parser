
import java.util.*;


public class Genes {
	String gene_id;
	String type;
	String strand; //use "+" and "-" or a boolean
	long gff_gene_start;
	long bed_gene_start;
	long gene_end;
	Chromosome chr; 
	Flanking up5k;
	Flanking dw5k;
	Set<Transcripts> transcripts = new HashSet<Transcripts>();
	List<Exon> exons = new ArrayList<Exon>();
	Set<Introns> introns = new HashSet<Introns>();
	Set<CDS> cds = new HashSet<CDS>();

	Genes(String id, String t, String s, long st, long e, Chromosome c ) {
		gene_id = id;
		type = t;
		strand = s;
		gff_gene_start = st;
		bed_gene_start = st - 1;
		gene_end = e;
		chr = c;
		
		if(strand.equalsIgnoreCase("-")){
			long gffstartup = gene_end;
			long gffendup = movecoord(chr, gene_end, 5000, true );
			long gffstartdown = movecoord(chr, gff_gene_start, 5000, false );
			long gffenddown = gff_gene_start;
			up5k = new Flanking (this, "Upstream", gffstartup, gffendup, "-", chr); 
			dw5k = new Flanking (this, "Downstream", gffstartdown, gffenddown, "-", chr);
		
		} else {
			long gffstartup = movecoord(chr, gff_gene_start, 5000, false );
			long gffendup = gff_gene_start; 
			long gffstartdown = gene_end;
			long gffenddown = movecoord(chr, gene_end, 5000, true );
			up5k = new Flanking (this, "Upstream", gffstartup, gffendup, "+", chr); 
			dw5k = new Flanking (this, "Downstream", gffstartdown, gffenddown, "+", chr); 
		} 
		
	}
	
	
	public long movecoord(Chromosome c, long start, long move, boolean add) {
		long returncoord = 0;
		if(add){
			long coord = start + move;
			if(coord > c.getChrsize()){
				returncoord = c.getChrsize();
			} else {
				returncoord = coord;
			}
			
		} else {
			long coord = start - move;
			if(coord < 1) {
				returncoord = 1;
			} else {
				returncoord = coord;
			}
		}
		return returncoord;
	}
	
	public void addTranscripts(Transcripts trans) {
		transcripts.add(trans);
	}
	
	public void addExons(Exon ex){
		exons.add(ex);
	}
	
	public void addIntrons(Introns in){
		introns.add(in);
	}
	
	public void addCDSs(CDS c){
		cds.add(c);
	}
	
	public boolean isproteincoding() {
		if(type.equalsIgnoreCase("protein_coding")){
			return true;
		} else {
			return false;
		}
	}
	
	public boolean ispseudogene() {
		if(type.equalsIgnoreCase("pseudogene")){
			return true;
		} else {
			return false;
		}
	}
	
	public boolean istransposon() {
		if(type.equalsIgnoreCase("transposable_element")){
			return true;
		} else {
			return false;
		}
	}

	public String getGene_id() {
		return gene_id;
	}

	public long getGff_gene_start() {
		return gff_gene_start;
	}

	public long get_gene_end() {
		return gene_end;
	}

	public Chromosome getChr() {
		return chr;
	}

	public Flanking getUp5k() {
		return up5k;
	}

	public Flanking getDw5k() {
		return dw5k;
	}
	
	public void setUp5k(Flanking f){
		up5k = f;		
	}
	
	public void setDw5k(Flanking f){
		dw5k = f;		
	}
	
}