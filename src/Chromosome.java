import java.util.ArrayList;
import java.util.List;

public class Chromosome {
		String chrname = "";
		int chr; //if chromosome is unknown = 0, if chromosome is mt = 1000 and if it is pt = 100
		long chrstart;
		long chrend;
		long chrsize;
		List<Genes> chrgenelist = new ArrayList<Genes>();
		List<Genes> chrpseudolist = new ArrayList<Genes>();
		List<Genes> chrtranspolist = new ArrayList<Genes>();
        
		Chromosome(String name, int nm, long start, long end) {
			chrname = name;
			chr = nm;
			chrstart = start;
			chrend = end;
			chrsize = end;
		}

		public List<Genes> getChrgenelist() {
			return chrgenelist;
		}
		
		public List<Genes> getpseudogenes() {
			return chrpseudolist;
		}
		
		public List<Genes> gettransposons() {
			return chrtranspolist;
		}

		public void setChrgenelist(List<Genes> chrgenelist) {
			this.chrgenelist = chrgenelist;
		}
		
		public long getChrsize(){
			return chrsize;
		}
		
		public boolean isNuclear(){
			if(chr > 0 && chr < 99) {
				return true;
			} else {
				return false;
			}
		}
		
		public boolean isUNKNOWN() {
			if(chr == 0){
				return true;
			} else {
				return false;
			}
		}
		
		public boolean isMitochondrial(){
			if (chr == 100) {
				return true;
			} else {
				return false;
			}
		}
		
		public boolean isPlastid(){
			if (chr == 1000) {
				return true;
			} else {
				return false;
			}
		}
		
    }