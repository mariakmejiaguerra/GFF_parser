

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.text.ParseException;
import java.util.Hashtable;

//TODO - Check format for gff file and give exceptions or options to manage it
public class ParseGFF {
	
	static String gff_file;
	Hashtable<String, Chromosome> chromhash = new Hashtable<String, Chromosome>();
	Hashtable<String, Genes> allgenes = new Hashtable<String, Genes>();
	Hashtable<String, Transcripts> allmrnas = new Hashtable<String, Transcripts>();
	Hashtable<String, Genes> proteinhash = new Hashtable<String, Genes>(); //protein_coding
	Hashtable<String, Genes> transpohash = new Hashtable<String, Genes>(); //transposable_element
	Hashtable<String, Genes> pseudohash = new Hashtable<String, Genes>(); //pseudogene
	Hashtable<String, Transcripts> mrnahash = new Hashtable<String, Transcripts>();
	Hashtable<String, Flanking> Up5Kb = new Hashtable<String, Flanking>();
	Hashtable<String, Flanking> Down5Kb = new Hashtable<String, Flanking>();
	
	public ParseGFF(String fileName) throws IOException, ParseException {
		int count = 0;
		
		BufferedReader br = new BufferedReader(new FileReader(fileName));
		String line;
		while ((line = br.readLine()) != null) {
			count++;
			String[] items = line.split("\t");

			if (items[2].contentEquals("chromosome")) {
				if (items[0].equalsIgnoreCase("UNKNOWN")) {
						int chrom = 0;
						int startchr = Integer.parseInt(items[3].trim());
						int endchr = Integer.parseInt(items[4].trim());
						Chromosome mychrom = new Chromosome (items[0].trim() , chrom , startchr , endchr );
						chromhash.put(items[0].trim(), mychrom);
				} else if (items[0].equalsIgnoreCase("Pt")) {
						int chrom = 1000;
						int startchr = Integer.parseInt(items[3].trim());
						int endchr = Integer.parseInt(items[4].trim());
						Chromosome mychrom = new Chromosome (items[0].trim() , chrom , startchr , endchr );
						chromhash.put(items[0].trim(), mychrom);
						
				} else if (items[0].equalsIgnoreCase("Mt")) {
						int chrom = 100;
						int startchr = Integer.parseInt(items[3].trim());
						int endchr = Integer.parseInt(items[4].trim());
						Chromosome mychrom = new Chromosome (items[0].trim() , chrom , startchr , endchr );
						chromhash.put(items[0].trim(), mychrom);
										
				} else {
						int chrom = Integer.parseInt(items[0].trim());;
						int startchr = Integer.parseInt(items[3].trim());
						int endchr = Integer.parseInt(items[4].trim());
						Chromosome mychrom = new Chromosome (items[0].trim() , chrom , startchr , endchr);
						chromhash.put(items[0].trim(), mychrom);
						
				}
				
			} else if (items[2].equalsIgnoreCase("gene")) {
				if (items[6].equalsIgnoreCase("-")) {
					int startgene = Integer.parseInt(items[3].trim()) ;
					int endgene = Integer.parseInt(items[4].trim());
					Chromosome genechr = chromhash.get(items[0].trim());
					String[] decoration = items[8].split(";");
					String[] geneid = decoration[0].split("=");
					String[] biotype = decoration[2].split("=");
					Genes mygene = new Genes(geneid[1].trim(), biotype[1].trim(), "-", startgene, endgene, genechr);
					allgenes.put(geneid[1].trim(), mygene);
					setHashCollection(genechr, mygene);
												
				} else {
					int startgene = Integer.parseInt(items[3].trim()) ;
					int endgene = Integer.parseInt(items[4].trim());
					Chromosome genechr = chromhash.get(items[0].trim());
					String[] decoration = items[8].split(";");
					String[] geneid = decoration[0].split("=");
					String[] biotype = decoration[2].split("=");
					Genes mygene = new Genes(geneid[1].trim(), biotype[1].trim(), "+", startgene, endgene, genechr);
					allgenes.put(geneid[1].trim(), mygene);
					setHashCollection(genechr, mygene);
																		
				}
				
			} else if(items[2].equalsIgnoreCase("mRNA")) {
				int startmRNA = Integer.parseInt(items[3].trim());
				int endmRNA = Integer.parseInt(items[4].trim());
				String mRNAchr = items[0].trim();
				String[] decoration = items[8].split(";");
				String[] transcriptid = decoration[0].split("=");
				String[] parentgene = decoration[1].split("=");
				Transcripts myRNA = new Transcripts(transcriptid[1].trim(), parentgene[1].trim(), items[6].trim(), startmRNA, endmRNA, mRNAchr);
				allgenes.get(parentgene[1].trim()).addTranscripts(myRNA);
				allmrnas.put(myRNA.transcript_id, myRNA);
				
					
			} else if (items[2].equalsIgnoreCase("exon")){
				int startexon = Integer.parseInt(items[3].trim());
				int endexon = Integer.parseInt(items[4].trim());
				String exonchr = items[0].trim();
				String[] decoration = items[8].split(";");
				String[] exonid = decoration[0].split("=");
				String[] parentmrna = decoration[1].split("=");
				Exon exon = new Exon(parentmrna[1].trim(), exonid[1].trim(),  items[6].trim(), startexon, endexon, exonchr);
				allmrnas.get(exonid[1].trim()).addExons(exon);
				allgenes.get(allmrnas.get(exonid[1].trim()).parent_gene_id).addExons(exon);
								
			} else if (items[2].equalsIgnoreCase("intron")){
				int startintron = Integer.parseInt(items[3].trim());
				int endintron = Integer.parseInt(items[4].trim());
				String intronchr = items[0].trim();
				String[] decoration = items[8].split(";");
				String[] intronid = decoration[0].split("=");
				String[] parentmrna = decoration[1].split("=");
				Introns intron = new Introns(parentmrna[1].trim(), intronid[1].trim(), items[6].trim(), startintron, endintron, intronchr);
				allmrnas.get(intronid[1].trim()).addIntrons(intron);
				allgenes.get(allmrnas.get(intronid[1].trim()).parent_gene_id).addIntrons(intron);
								
			} else if (items[2].equalsIgnoreCase("CDS")) {
				int startCDS = Integer.parseInt(items[3].trim());
				int endCDS = Integer.parseInt(items[4].trim());
				String CDSchr = items[0].trim();
				String[] decoration = items[8].split(";");
				String[] CDSid = decoration[0].split("=");
				String[] parentmrna = decoration[1].split("=");
				CDS codingseq = new CDS( parentmrna[1].trim(), CDSid[1].trim(), items[6].trim(), startCDS, endCDS, CDSchr);
				allmrnas.get(CDSid[1].trim()).addCDS(codingseq);
				allgenes.get(allmrnas.get(CDSid[1].trim()).parent_gene_id).addCDSs(codingseq);
			}
		}
		
		br.close();
	}	
		
	public void setHashCollection(Chromosome c, Genes mg) {
		Chromosome genechr = c;
		Genes mygene = mg;
		if (genechr.isNuclear()) {
			if (mygene.isproteincoding()){
				genechr.getChrgenelist().add(mygene);
				proteinhash.put(mygene.getGene_id(), mygene);	
				Up5Kb.put(mygene.getGene_id(), mygene.getUp5k()); 
				Down5Kb.put(mygene.getGene_id(), mygene.getDw5k()); 
			} else if (mygene.ispseudogene()) {
				genechr.getpseudogenes().add(mygene);
				pseudohash.put(mygene.getGene_id(), mygene);
			} else if (mygene.istransposon()) {
				genechr.gettransposons().add(mygene);
				transpohash.put(mygene.getGene_id(), mygene);							
			}
		}
		
		if (genechr.isPlastid()) {
			if (mygene.isproteincoding()){
				chromhash.get("Pt").getChrgenelist().add(mygene);
				proteinhash.put(mygene.getGene_id(), mygene);
			} else if (mygene.ispseudogene()) {
				genechr.getpseudogenes().add(mygene);
				pseudohash.put(mygene.getGene_id(), mygene);
			} else if (mygene.istransposon()) {
				genechr.gettransposons().add(mygene);
				transpohash.put(mygene.getGene_id(), mygene);
			}
		} else if (genechr.isMitochondrial()) {
			if(mygene.isproteincoding()){
				chromhash.get("Mt").getChrgenelist().add(mygene);
				proteinhash.put(mygene.getGene_id(), mygene);
			} else if (mygene.ispseudogene()) {
				genechr.getpseudogenes().add(mygene);
				pseudohash.put(mygene.getGene_id(), mygene);
			} else if (mygene.istransposon()) {
				genechr.gettransposons().add(mygene);
				transpohash.put(mygene.getGene_id(), mygene);
			}
		} else if (genechr.isUNKNOWN()) {
			if(mygene.isproteincoding()) { 
				chromhash.get("UNKNOWN").getChrgenelist().add(mygene);
				proteinhash.put(mygene.getGene_id(), mygene);
			} else if (mygene.ispseudogene()) {
				genechr.getpseudogenes().add(mygene);
				pseudohash.put(mygene.getGene_id(), mygene);
			} else if (mygene.istransposon()) {
				genechr.gettransposons().add(mygene);
				transpohash.put(mygene.getGene_id(), mygene);
			}	
		}
		
	}
	
	@SuppressWarnings("unchecked")
	public Hashtable gethash(String evaluate){
		String asked = evaluate;
		if(asked.equalsIgnoreCase("chromosomes")) {
			return chromhash;
		} else if (asked.equalsIgnoreCase("protein")) {
			return proteinhash;
		} else if (asked.equalsIgnoreCase("transposons")) {
			return transpohash;
		} else if (asked.equalsIgnoreCase("pseudogenes")) {
			return pseudohash;
		} else if (asked.equalsIgnoreCase("transcripts")) {
			return mrnahash;
		} else if (asked.equalsIgnoreCase("Up5Kb")) {
			return Up5Kb;
		} else if (asked.equalsIgnoreCase("Down5Kb")) {
			return Down5Kb;
		}		
		return allgenes;
	}
	
	public static void main(String[] args) throws Exception {
		if (args.length != 1) {
			System.out.println("Usage: ParseGFF option {gff_file}");
			return;
		}
		gff_file = args[0];
	    System.out.println("Parse Report Version 1 (February 2011)");
		ParseGFF obj = new ParseGFF(gff_file);
		
		@SuppressWarnings("unused")
		Output out_to_beds = new Output(obj);
		out_to_beds.make_genesbedfile("Zm5b60_fromparsegff.bed");
		//out_to_beds.make_bedfile_pergene("GRMZM2G422750", 3000);
	}
}