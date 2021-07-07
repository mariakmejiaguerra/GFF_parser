


public class Exon {
	String gffexon_id;
	String exon_id;
	String parent_transcript_id;
	String strand; //use "+" and "-" or a boolean
	String exonchr;
	long gff_exon_start;
	long bed_exon_start;
	long exon_end;

	public Exon(String id, String pid, String str, long st, long ed, String c) {
		gffexon_id = id;
		parent_transcript_id = pid;
		strand = str;
		gff_exon_start = st;
		bed_exon_start = st - 1;
		exon_end = ed;
		exonchr = c;
	}
	
	public void setExonID(String id) {
		exon_id = id;
	}
}