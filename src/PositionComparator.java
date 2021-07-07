

import java.util.Comparator;


public class PositionComparator implements Comparator<Object> {

boolean negativestrand;	
	public PositionComparator(boolean strand) {
		negativestrand = strand;
	}
	
	public int compare(Object o1, Object o2) {
		if (o1 instanceof Exon && o2 instanceof Exon) {
			return compareExons((Exon)o1, (Exon)o2);
		} else if (o1 instanceof CDS && o2 instanceof CDS) {
			return compareCDS((CDS)o1, (CDS)o2);
		} else if (o1 instanceof Introns && o2 instanceof Introns) {
			return compareIntrons((Introns)o1, (Introns)o2);
		}
		return 0;
	}
	
	public int compareExons(Exon e1, Exon e2) {
		if (negativestrand){
			if( e1.exon_end > e2.exon_end ) {
				return -1;
			} else if ( e1.exon_end < e2.exon_end ) {
				return 1;
			}   
			return 0;
		} else {
			if( e1.gff_exon_start > e2.gff_exon_start ) {
				return 1;
			} else if ( e1.gff_exon_start < e2.gff_exon_start ) {
				return -1;
			}   
			return 0;			
		}
	}
	
	public int compareIntrons (Introns i1, Introns i2) {
		if (negativestrand){
			if( i1.intron_end > i2.intron_end ) {
				return -1;
			} else if ( i1.intron_end < i2.intron_end ) {
				return 1;
			}   
		return 0;
		} else {
			if( i1.gff_intron_start > i2.gff_intron_start ) {
				return 1;
			} else if( i1.gff_intron_start < i2.gff_intron_start ) {
				return -1;
			}   
		return 0;
		}
	}
	
	public int compareCDS (CDS cds1, CDS cds2) {
		if(negativestrand){
			if (cds1.CDS_end > cds2.CDS_end) {
				return -1;
			} else if (cds1.CDS_end < cds2.CDS_end) {
				return 1;
			} return 0;
		} else {
			if (cds1.gff_CDS_start > cds2.gff_CDS_start) {
				return 1;
			} else if (cds1.gff_CDS_start < cds2.gff_CDS_start) {
				return -1;
			} return 0;
		}
		
		
	}

	

}
