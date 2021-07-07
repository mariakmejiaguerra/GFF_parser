import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Collections;
import java.util.Hashtable;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.Set;

public class Output {
	
	Hashtable<String, Chromosome> chromosomes = new Hashtable<String, Chromosome>();
	Hashtable<String, Genes> allgenes = new Hashtable<String, Genes>();
	Hashtable<String, Genes> proteingenes = new Hashtable<String, Genes>(); //protein_coding
	Hashtable<String, Genes> transpogenes = new Hashtable<String, Genes>(); //transposable_element
	Hashtable<String, Genes> pseudogenes = new Hashtable<String, Genes>(); //pseudogene
	Hashtable<String, Transcripts> mrnas = new Hashtable<String, Transcripts>();
	Hashtable<String, Flanking> FlankingUp5Kb = new Hashtable<String, Flanking>();
	Hashtable<String, Flanking> FlankingDown5Kb = new Hashtable<String, Flanking>();
	
	@SuppressWarnings("unchecked")
	public Output(ParseGFF parsegff) {
		chromosomes = parsegff.gethash("chromosomes");
		proteingenes = parsegff.gethash("protein");
		transpogenes = parsegff.gethash("transposons");
		pseudogenes = parsegff.gethash("pseudogenes");
		mrnas = parsegff.gethash("transcripts");
		allgenes = parsegff.gethash("all");
		//FlankingUp5Kb = parsegff.gethash("Up5Kb");
		//FlankingDown5Kb = parsegff.gethash("Down5Kb");	
	}
	
	public void make_genesbedfile(String file_name) throws IOException{
		FileWriter outFile = new FileWriter(file_name);
		PrintWriter out = new PrintWriter(outFile);		
		StringBuffer bedfile = new StringBuffer();
		Set<String> idn = allgenes.keySet();
		Iterator<String> itr1 = idn.iterator();
		while ( itr1.hasNext() ) {
			String hashkey = itr1.next();
			Genes eachgene = allgenes.get(hashkey);
			Chromosome chromnow = eachgene.getChr();
			
			bedfile.append(chromnow.chrname+"\t"+eachgene.bed_gene_start+"\t"+eachgene.gene_end+"\t"+eachgene.gene_id+"\t"+"."+"\t"+eachgene.strand+"\t"+eachgene.type+"\n");
			bedfile.append(chromnow.chrname+"\t"+eachgene.getUp5k().bed_flank_start+"\t"+eachgene.getUp5k().flank_end+
					"\t"+"Up5kb_"+eachgene.gene_id+"\t"+"."+"\t"+eachgene.strand+"\t"+"flanking_region"+"\n");
			bedfile.append(chromnow.chrname+"\t"+eachgene.getDw5k().bed_flank_start+"\t"+eachgene.getDw5k().flank_end+
					"\t"+"Down5kb_"+eachgene.gene_id+"\t"+"."+"\t"+eachgene.strand+"\t"+"flanking_region"+"\n");
			Set<Transcripts> mrnas = eachgene.transcripts;
			for (Transcripts eachmrnas : mrnas) {
				bedfile.append(chromnow.chrname+"\t"+eachmrnas.bed_mrna_start+"\t"+eachmrnas.mrna_end+"\t"+eachmrnas.transcript_id+"\t"+"."+"\t"+eachmrnas.strand+"\t"+"transcript"+"\n");
				List<Exon> exons = new LinkedList<Exon>( eachmrnas.getExons() );
				List<CDS> cdss = new LinkedList<CDS>(eachmrnas.getCDS());
				List<Introns> introns = new LinkedList<Introns>(eachmrnas.getIntrons());
				Collections.sort( exons, new PositionComparator(eachmrnas.strand.equalsIgnoreCase("-")));
				Collections.sort( cdss, new PositionComparator(eachmrnas.strand.equalsIgnoreCase("-")));
				Collections.sort( introns, new PositionComparator(eachmrnas.strand.equalsIgnoreCase("-")));
				int e = 0;
				for (Exon eachexon : exons){
					bedfile.append(chromnow.chrname+"\t"+ eachexon.bed_exon_start+"\t"+ eachexon.exon_end+"\t"+eachexon.gffexon_id+"\t"+"."+"\t"+eachexon.strand+"\t"+"exon"+e+"\n");
					e++;
				}
				int i = 0;
				for (Introns eachintron : introns) {
					bedfile.append(chromnow.chrname+"\t"+eachintron.bed_intron_start+"\t"+eachintron.intron_end+"\t"+eachintron.intron_id+"\t"+"."+"\t"+eachintron.strand+"\t"+"intron"+i+"\n");
					i++;
				}
				int c = 0;
				for (CDS eachcds : cdss) {
					bedfile.append(chromnow.chrname+"\t"+eachcds.bed_CDS_start+"\t"+eachcds.CDS_end+"\t"+eachcds.CDS_id+"\t"+"."+"\t"+eachcds.strand+"\t"+"CDS"+c+"\n");
					c++;
				}
			}
		}
		out.println(bedfile);
		out.close();
	}
	
	public void make_bedfile_pergene(String gene_name, int flanking) throws IOException {
		String hashkey = gene_name.trim();
		FileWriter outFile = new FileWriter(gene_name+".bed");
		PrintWriter out = new PrintWriter(outFile);
		StringBuffer bedfile = new StringBuffer();
		Genes eachgene = allgenes.get(hashkey);
		Chromosome chromnow = eachgene.getChr();
		long upstreamstart;
		long upstreamend;
		long downstreamstart;
		long downstreamend;
		if(eachgene.strand.equals("-")){
			upstreamstart = eachgene.gene_end;
			upstreamend = eachgene.movecoord(chromnow, eachgene.gene_end, flanking, true);
			downstreamstart = eachgene.movecoord(chromnow, eachgene.bed_gene_start, flanking, false);
			downstreamend = eachgene.bed_gene_start;
		} else {
			upstreamend  = eachgene.bed_gene_start; 
			downstreamstart = eachgene.gene_end;
			upstreamstart = eachgene.movecoord(chromnow, eachgene.bed_gene_start, flanking, false );
			downstreamend = eachgene.movecoord(chromnow, eachgene.gene_end, flanking, true );
		}
		
		bedfile.append(chromnow.chrname+"\t"+eachgene.bed_gene_start+"\t"+eachgene.gene_end+"\t"+eachgene.gene_id+"\t"+"."+"\t"+eachgene.strand+"\t"+eachgene.type+"\n");
		bedfile.append(chromnow.chrname+"\t"+upstreamstart+"\t"+upstreamend+"\t"+"5'flanking_"+eachgene.gene_id+"\t"+"."+"\t"+eachgene.strand+"\t"+"flanking_upstream"+"\n");
		bedfile.append(chromnow.chrname+"\t"+downstreamstart+"\t"+downstreamend+"\t"+"3'flanking_"+eachgene.gene_id+"\t"+"."+"\t"+eachgene.strand+"\t"+"flanking_downstream"+"\n");
		eachgene.movecoord(chromnow, eachgene.bed_gene_start, 3000, true);
		Set<Transcripts> mrnas = eachgene.transcripts;
		for (Transcripts eachmrnas : mrnas) {
			bedfile.append(chromnow.chrname+"\t"+eachmrnas.bed_mrna_start+"\t"+eachmrnas.mrna_end+"\t"+eachmrnas.transcript_id+"\t"+"."+"\t"+eachmrnas.strand+"\t"+"transcript"+"\n");
			List<Exon> exons = new LinkedList<Exon>( eachmrnas.getExons() );
			List<CDS> cdss = new LinkedList<CDS>(eachmrnas.getCDS());
			List<Introns> introns = new LinkedList<Introns>(eachmrnas.getIntrons());
			Collections.sort( exons, new PositionComparator(eachmrnas.strand.equalsIgnoreCase("-")));
			Collections.sort( cdss, new PositionComparator(eachmrnas.strand.equalsIgnoreCase("-")));
			Collections.sort( introns, new PositionComparator(eachmrnas.strand.equalsIgnoreCase("-")));
			for (Exon eachexon : exons){			
				bedfile.append(chromnow.chrname+"\t"+ eachexon.bed_exon_start+"\t"+ eachexon.exon_end+"\t"+eachexon.gffexon_id+"\t"+"."+"\t"+eachexon.strand+"\t"+"exon"+"\n");
			} 
			for (Introns eachintron : introns) {
				bedfile.append(chromnow.chrname+"\t"+eachintron.bed_intron_start+"\t"+eachintron.intron_end+"\t"+eachintron.intron_id+"\t"+"."+"\t"+eachintron.strand+"\t"+"intron"+"\n");
			}
			for (CDS eachcds : cdss) {
				bedfile.append(chromnow.chrname+"\t"+eachcds.bed_CDS_start+"\t"+eachcds.CDS_end+"\t"+eachcds.CDS_id+"\t"+"."+"\t"+eachcds.strand+"\t"+"CDS"+"\n");
			}
		}
		out.println(bedfile);
		out.close();		
	}
	
	public void make_bedfile_pergene(String gene_name) throws IOException {
		String hashkey = gene_name.trim();
		FileWriter outFile = new FileWriter(gene_name);
		PrintWriter out = new PrintWriter(outFile);
		StringBuffer bedfile = new StringBuffer();
		Genes eachgene = allgenes.get(hashkey);
		Chromosome chromnow = eachgene.getChr();
			
		bedfile.append(chromnow.chrname+"\t"+eachgene.bed_gene_start+"\t"+eachgene.gene_end+"\t"+eachgene.gene_id+"\t"+"."+"\t"+eachgene.strand+"\t"+eachgene.type+"\n");
		bedfile.append(chromnow.chrname+"\t"+eachgene.getUp5k().bed_flank_start+"\t"+eachgene.getUp5k().flank_end+
				"\t"+"Up5kb_"+eachgene.gene_id+"\t"+"."+"\t"+eachgene.strand+"\t"+"flanking_region"+"\n");
		bedfile.append(chromnow.chrname+"\t"+eachgene.getDw5k().bed_flank_start+"\t"+eachgene.getDw5k().flank_end+
				"\t"+"Down5kb_"+eachgene.gene_id+"\t"+"."+"\t"+eachgene.strand+"\t"+"flanking_region"+"\n");
		Set<Transcripts> mrnas = eachgene.transcripts;
		for (Transcripts eachmrnas : mrnas) {
			bedfile.append(chromnow.chrname+"\t"+eachmrnas.bed_mrna_start+"\t"+eachmrnas.mrna_end+"\t"+eachmrnas.transcript_id+"\t"+"."+"\t"+eachmrnas.strand+"\t"+"transcript"+"\n");
			List<Exon> exons = new LinkedList<Exon>( eachmrnas.getExons() );
			List<CDS> cdss = new LinkedList<CDS>(eachmrnas.getCDS());
			List<Introns> introns = new LinkedList<Introns>(eachmrnas.getIntrons());
			Collections.sort( exons, new PositionComparator(eachmrnas.strand.equalsIgnoreCase("-")));
			Collections.sort( cdss, new PositionComparator(eachmrnas.strand.equalsIgnoreCase("-")));
			Collections.sort( introns, new PositionComparator(eachmrnas.strand.equalsIgnoreCase("-")));
			for (Exon eachexon : exons){			
				bedfile.append(chromnow.chrname+"\t"+ eachexon.bed_exon_start+"\t"+ eachexon.exon_end+"\t"+eachexon.gffexon_id+"\t"+"."+"\t"+eachexon.strand+"\t"+"exon"+"\n");
			} 
			for (Introns eachintron : introns) {
				bedfile.append(chromnow.chrname+"\t"+eachintron.bed_intron_start+"\t"+eachintron.intron_end+"\t"+eachintron.intron_id+"\t"+"."+"\t"+eachintron.strand+"\t"+"intron"+"\n");
			}
			for (CDS eachcds : cdss) {
				bedfile.append(chromnow.chrname+"\t"+eachcds.bed_CDS_start+"\t"+eachcds.CDS_end+"\t"+eachcds.CDS_id+"\t"+"."+"\t"+eachcds.strand+"\t"+"CDS"+"\n");
			}
		}
		out.println(bedfile);
		out.close();		
	}
	
	public void make_genomefile (String file_name) throws IOException {
		FileWriter outFile = new FileWriter(file_name);
		PrintWriter out = new PrintWriter(outFile);
		StringBuffer genomefile = new StringBuffer();
		for( String key:  chromosomes.keySet() ) {
			Chromosome chromnow = chromosomes.get(key);
			long endbound = chromnow.chrend;
			genomefile.append(chromnow.chrname+"\t"+endbound+"\n");
		}
		out.println(genomefile);
		out.close();
	}
	
	/**Set<Transcripts> mrnas = eachgene.transcripts;
	for (Transcripts eachmrnas : mrnas) {
		long tstart = eachmrnas.gff_mrna_start;
		long tend = eachmrnas.mrna_end;
		bedfile.append(chromnow.chrname+"\t"+tstart+"\t"+tend+"\t"+eachgene.gene_id+"\t"+"."+"\t"+eachgene.strand+"\n");
		List<Exon> exons = new LinkedList<Exon>( eachmrnas.getExons() );
		Collections.sort( exons, new PositionComparator(eachmrnas.strand.equalsIgnoreCase("-")) );					
	} **/
}
