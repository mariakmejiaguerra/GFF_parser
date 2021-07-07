


public class CDS {
	String CDS_id;
	String parent_transcript_id;
	String strand; //use "+" and "-" or a boolean
	String CDSchr;
	long gff_CDS_start;
	long bed_CDS_start;
	long CDS_end;
	
	public CDS(String id, String pid, String str, long st, long ed, String c) {
		CDS_id = id;
		parent_transcript_id = pid;
		strand = str;
		gff_CDS_start = st;
		bed_CDS_start = st - 1;
		CDS_end = ed;
		CDSchr = c;
	}
}