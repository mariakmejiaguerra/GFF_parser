


public class Flanking {
	Genes gene_id;
	String type;
	String strand;
	long gff_flank_start;
	long bed_flank_start;
	long flank_end;
	Chromosome chr;
	
	
	Flanking (Genes id, String tp, long start, long end, String s, Chromosome chr){
		gene_id = id;
		type = tp;
		gff_flank_start = start;
		if (start == 1 || start > 1) {
			bed_flank_start = start - 1;
		}		
		flank_end = end;
		strand = s;
	}

	public long getStart() {
		return gff_flank_start;
	}

	public long getEnd() {
		return flank_end;
	}
	
	public String getStrand() {
		return strand;
	}
	
	public boolean getType(){
		if (type.equalsIgnoreCase("Upstream")){
			return true;
		} else {
			return false;
		}
	}

}
