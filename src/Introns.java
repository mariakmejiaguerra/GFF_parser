


public class Introns {
	String intron_id;
	String parent_transcript_id;
	String strand; //use "+" and "-" or a boolean
	String intronchr;
	long gff_intron_start;
	long bed_intron_start;
	long intron_end;
	
	public Introns(String id, String pid, String str, long st, long ed, String c) {
		intron_id = id;
		parent_transcript_id = pid;
		strand = str;
		gff_intron_start = st;
		bed_intron_start = st -1 ;
		intron_end = ed;
		intronchr = c;
	}
	
	
	

}