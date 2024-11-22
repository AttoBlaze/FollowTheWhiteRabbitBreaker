import java.io.BufferedReader;
import java.io.FileReader;
import java.math.BigInteger;
import java.security.MessageDigest;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.function.Predicate;

public class TrustPilotMD5 {
	/*
	 * This code is literally unreadable but its fast asf so have at you
	 * 
	 * First, words are found by checking if they contain more of each letter than the anagram has
	 * Then, all the different ways you can achieve the length of the anagram by adding specific words lengths are found (coinchange)
	 * Then, all distinct permutations of sentences (aka. sentences impartial to word order) are tested to see if they are anagrams
	 * Checking if something is an anagram is done by assigning a prime to each letter and getting the prime factorization product of all words and checking if it is the same as the anagram factorization
	 * Then, all possible permutations (aka all possible anagrams) are made by permutating every distinct permutation.
	 * Lastly, all anagrams are checked to see if their hash is the same as the specified hash.
	 * This is all run on multiple threads to optimize.
	 * 
	 * */
	
	
	private static String[] wordList;
	private static String anagramString = "poultry outwits ants";
	private static int[] anagramCharList = anagramToCharacterCount(anagramString.replaceAll("\\s+",""));
	private static HashMap<Integer, String[]> wordByLength = new HashMap<Integer,String[]>();
	private static String[][] wordByLengthArray;
	private static long[][] wordFactorizationByLengthArray;
	private static int[] wordLengths;
	private static HashSet<String> anagramsFound;
	private static final int threads = 16;
	private static final int statusInterval = 1000;
	private static long anagramFactorization;
	private static int[] primes;
	private static String messageHash = "";
	private static ArrayList<Long> tests = new ArrayList<Long>();
	private static long startTime;
	
	private static String Easy = 	"e4820b45d2277f3844eac66c903e84be",
						  Medium = 	"23170acc097c24edb98fc5488ab033fe",
						  Hard = 	"665e5bcb0c20062fe8abaaf4628bb154";
	
	//converts a string to a list of character frequencies
	public static int[] anagramToCharacterCount(String anagram) {
		int[] charCount = new int[26];
		
		//get count of all characters
		for(int i=0;i<anagram.length();i++) {
			char c = anagram.charAt(i);
			charCount[c-'a']++;
		}
		return charCount;
	}
	
	//checks if a number is prime
	public static boolean isPrime(int n) {
		for (int i=2;i<Math.sqrt(n);i++) {
			if (n%i==0) return false;
		}
		return true;
	}
	
	//gets prime factorization according to how many of each prime number exists
	private static long GetFactorization(int[] primeCounts) {
		long result = 1;
		for (int i=0;i<primeCounts.length;i++) for (int j=0;j<primeCounts[i];j++) result*=primes[i];
		return result;
	}
	
	public TrustPilotMD5() {
		//get primes
		ArrayList<Integer> prime = new ArrayList<Integer>();
		for(int i=2 ; prime.size()<anagramCharList.length ; i++) if(isPrime(i)) prime.add(i);
		
		//sort primes by frequency in anagram to prevent overflow when factorizing (for matches)
		class primeFrequencySorter{
			public int index;
			public int value;
			public primeFrequencySorter(int a, int b) {index=a;value=b;}
		}
		ArrayList<primeFrequencySorter> sorter = new ArrayList<primeFrequencySorter>();
		for (int i=0;i<prime.size();i++) sorter.add(new primeFrequencySorter(i,anagramCharList[i]));
		class thisCompare implements Comparator<primeFrequencySorter>{
			public int compare(primeFrequencySorter i, primeFrequencySorter j) {
				return j.value-i.value;
		}}
		sorter.sort(new thisCompare());
		primes = new int[prime.size()];
		for (int i=0;i<prime.size();i++) primes[sorter.get(i).index]=prime.get(i);
		
		//get anagram factorization
		anagramFactorization = GetFactorization(anagramCharList);
		  
		//read words file to find possible anagram words
		startTime = System.nanoTime();
		ReadWordFile();
		System.out.println("Added " + wordList.length + " possible words in "+(float)((System.nanoTime()-startTime)/1000)/1000+"ms");
		
		messageHash = Easy;
		System.out.println("\nEasy hash (3 words):");
		CompleteTest(x -> x==3);
		
		messageHash = Medium;
		System.out.println("\nMedium hash (3 words):");
		CompleteTest(x -> x==3);

		messageHash = Hard;
		System.out.println("\nHard hash (4 words):");
		CompleteTest(x -> x==4);
	}
	
	private static void CompleteTest(Predicate<Integer> countConstriction) {
		//set stuff up
        anagramsFound = new HashSet<String>();
        tests.clear();
        startTime = System.nanoTime();
        
        //get word length combinations
        ArrayList<int[]> wordCounts = CoinChange(anagramString.replaceAll("\\s+","").length(),wordLengths,countConstriction);
        
        //get all anagrams impartial to word order
        GetAllAnagrams(wordCounts,wordLengths,anagramsFound,threads);
        int matches = anagramsFound.size();
        System.out.print(GetSum(tests)+" combinations checked in "+
				(float)((System.nanoTime()-startTime)/1000000)/1000
				+"s, of which "
				+matches
				+" matches were identified, translating to ");

        //get all permutations of anagrams (all unique anagrams) (anagrams partial to word order)
        GetPermutationsOfMatch(anagramsFound);
        int anagrams = anagramsFound.size();
        System.out.print(anagrams+" unique anagrams, of which ");
        
        //get all hash matches 
		GetHashMatches(anagramsFound);
		String hashMatches = "";
		for (String i:anagramsFound) hashMatches+=i+", ";
		if (anagramsFound.size()>0) hashMatches = hashMatches.substring(0,hashMatches.length()-2);
		System.out.println(+anagramsFound.size()+" match the hash (the whole process taking "+(float)((System.nanoTime()-startTime)/1000000)/1000+"s):\n"
								+"[ "+hashMatches+" ]");
	}
	
	//prints the number of combinations (amount of unique anagrams impartial to word order) checked, time, and matches (anagrams) found.
	private static void PrintTests(int matches, int totalVariations) {
		System.out.println(GetSum(tests)/1000000 + " million combinations checked in "+(float)((System.nanoTime()-startTime)/1000000)/1000+"s ("+matches+" matches found) ("+tests.size()+"/"+totalVariations+" word length combinations checked)");
	}
	
	//gets the sum of a list of longs (normally the amount of tests performed)
	private static long GetSum(ArrayList<Long> longs) {
		long result = 0;
		for (long l:longs) result+=l;
		return result;
	}
	
	//gets the MD5 hash of a string
	private static String GetHash(String s) {
		String result = "";
		try {
			MessageDigest md = MessageDigest.getInstance("MD5");
			md.update(s.getBytes());
			result = new BigInteger(1,md.digest()).toString(16);
		} 
		catch (Exception e) {e.printStackTrace();}
		
		return result;
	}
	
	//filters out all strings whose hash is not the same as the message hash
	private static void GetHashMatches(HashSet<String> anagrams) {
		String[] Anagrams = anagrams.toArray(new String[anagrams.size()]);
		anagrams.clear();
		
		for(String s:Anagrams) {
			if(GetHash(s).equals(messageHash)) anagrams.add(s);
		}
	}
	
	//gets a list of words the message must be comprised of after which a list of preporatory operations are performed
	public static void ReadWordFile() {
		//read file
		FileReader fReader;
		try {
			fReader = new FileReader("wordlist.txt");
			BufferedReader bReader = new BufferedReader(fReader);
			String word;
			var words = new HashSet<String>();
			
			//get possible words
			while((word = bReader.readLine()) != null) {
				//check if the word contains more of a letter than is possible in the anagram
				if (LettersAvailable(word)) {
					//handle single letter edge cases 
					if (word.length()>1 || word.matches("[aoi]")) words.add(word);
				}
			}
			bReader.close();
			fReader.close();
			
			//convert result to array
			int i=0;
			wordList = new String[words.size()];
			for(String ele:words) wordList[i++]=ele;
			
			//get hashmap of words by length (eg. 5 -> all 5 letter words)
			var wordsByLength = new HashMap<Integer,ArrayList<String>>();
			for(String ele:wordList) {
				if (wordsByLength.containsKey(ele.length())) { //add word to list of same letter count words
					var elem = wordsByLength.get(ele.length());
					elem.add(ele);
					wordsByLength.put(ele.length(),elem);
				}
				else wordsByLength.put(ele.length(),new ArrayList<String>()); //create new list if it does not exist for that letter count
			}
			
			//reformat hashmap to use arrays instead of arraylist and additionally create an array of word lengths
			int max = 0; i=0;
			wordLengths = new int[wordsByLength.keySet().size()]; 
			for(int key:wordsByLength.keySet()) {
				//convert arraylist to list
				var list = wordsByLength.get(key);
				String[] List = new String[list.size()];
				for(int j=0;j<list.size();j++) {
					List[j]=list.get(j);
				}
				
				//put
				wordByLength.put(key,List);
				wordLengths[i++] = key;
				max = Math.max(max,key);
			}
			
			//turn hashmap into array and get all word factorizations
			wordByLengthArray = new String[max+1][];
			wordFactorizationByLengthArray = new long[max+1][];
			for(int j=1;j<=max;j++) {
				wordByLengthArray[j] = wordByLength.get(j);									//get word list
				wordFactorizationByLengthArray[j] = new long[wordByLengthArray[j].length];	//get factorizations
				for (int k=0;k<wordByLengthArray[j].length;k++) {
					wordFactorizationByLengthArray[j][k] = GetFactorization(anagramToCharacterCount(wordByLengthArray[j][k]));
				}
			}
			
			//sort word lengths array
			Arrays.sort(wordLengths);
			for(int j=0;j<wordLengths.length/2;j++) {
				Swap(wordLengths,j,wordLengths.length-1-j);
			}
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
	
	
	
	
	
	
	
	
	
	
	//checks if the given string contains more of any one letter than appears in the anagram
	private static int[] charList;
	private static boolean LettersAvailable(String word) {
		charList = new int[26];
		for (int i=0;i<word.length();i++) {
			//get char index in alphabet
			int charIndex = word.charAt(i)-'a';
			
			//if char is not in alphabet
			if (charIndex<0 || charIndex>=anagramCharList.length ||

			//or if the count of a letter is greather than the anagram letter count then the word cannot be part of the anagram.
				++charList[charIndex]>anagramCharList[charIndex]) return false;
		}
		return true;
	}
	
	
	/* Coin change algorithm:
	 * |	The idea is that for any given list of coins and a sum (using only positive integers), how many unique ways can 
	 * |	you add together coins to get the sum?
	 * |	fx. for coins = [3,2,1] and sum = 5:
	 * |		5 = 3+2
	 * |		5 = 3+1+1
	 * |		5 = 2+2+1
	 * |		5 = 2+1+1+1
	 * |		5 = 1+1+1+1+1
	 * |	eg. anwser = 5.
	 * 	
	 * This is an original dynamic solution to the problem which i implemented in a different project.
	 * Here, it is used to find all the different ways you can add together words of given lengths to get a final
	 * sentence which contains n letters.
	 * So fx. for the previous example, think "how many different ways can you make a sentence consisting of 5 letters, using 
	 * only 3, 2, and 1 letter words?".
	 * Specifically, think of how many of each word length you would need to get the desired result instead of thinking about the actual words themselves 
	 * (eg. one 3 letter word and one 2 letter word can together make a 5 letter sentence (result here would be [1,1,0])).
	 * We additionally apply a constriction here to the total number of words used (Fx the final sentence must consist of 3 words (eg. x -> x==3))
	 * In the end this means that we can use this algorithm to get a list of how many words of each length we we will need to make a final sentence which
	 *  1: consist of an appropriate number of words
	 *  2: contain the same amout of letters as our message anagram.
	 * This effectively means that we can skip any process related to checking the length of words and can jump straight to creating sentences and finding anagrams.
	 * 
	 * Final note: 
	 * The comments in this method are copied over and as such are contextual to coin change. Remember that the ideas discussed still apply, even though the 
	 * subject of this method's comments is that of coins.
	 */
	private static ArrayList<int[]> CoinChange(int amount,int[] coins,Predicate<Integer> countConstriction) {
		//initialize stuff
        int currentCoinIndex = 0,
            currentCoinTotal = 0;
        int[] amountOfCoins = new int[coins.length];
        var variations = new ArrayList<int[]>(); 
        
        while(true) {
            //if we have added all the coins we can then check for a complete sum
            if (currentCoinIndex>=coins.length) {
                //if we have the correct total then we have completed a variation
                if (currentCoinTotal==amount) {
                	//get the total amount of coins used in this variation
                	int coinCount = 0; for (int i:amountOfCoins) coinCount+=i;
                	
                	//check if count constrictions apply, and if so, add the the list of variationss
                	if (countConstriction.test(coinCount)) {
                		variations.add(amountOfCoins.clone());
                	}
                }

                //go backwards to find the smallest coin we can remove to create a new variation
                do if (currentCoinIndex<=0) {//if the smallest coin we can remove is bigger than the biggest coin, then we have checked all variations. 
            		//return matches
            		return variations;  
                }
                while (amountOfCoins[--currentCoinIndex]==0);

                //remove one of said coins
                amountOfCoins[currentCoinIndex]--;
                currentCoinTotal-=coins[currentCoinIndex];
            }
            else {
            	//increase coin count
                int coinAmount = (int)(amount-currentCoinTotal)/coins[currentCoinIndex];
                amountOfCoins[currentCoinIndex] = coinAmount;
                currentCoinTotal+= coins[currentCoinIndex]*coinAmount;
            }
            
            //find the next coin we can add as the biggest coin smaller than the last removed coin
            do if (++currentCoinIndex>=coins.length) break;
            while (amount-currentCoinTotal<coins[currentCoinIndex]);
        }
    }
	
	
	
	
	
	
	//gets all unique anagrams impartial of word order (aka. all unique combinations of words which result in an anagram)
	private static void GetAllAnagrams(ArrayList<int[]> wordCounts, int[] wordLengths, HashSet<String> matches, int threadCount) {
		ThreadGroup workers = new ThreadGroup("anagramFinders");
		
		//create status writer (writes running test results)
        var statusWriter = new Runnable() {
        	public void run() {
        		int i=0;
        		long time = System.currentTimeMillis();
    			while (workers.activeCount()>0) {	
        			//print status
        			PrintTests(matches.size(),wordCounts.size());
        			
        			//sleep thread
    				try {Thread.sleep((++i)*statusInterval-(System.currentTimeMillis()-time));} 
    				catch (InterruptedException e) {}
        		}
        	}
        };
        
        //split up into seperate tasks (by giving out variations from coin change to different threads and checking each variation for anagrams)
		for (int i=0;i<threadCount;i++) {
			int k = i;
			//create new thread with task
    		Runnable GetSinglePermutatonAnagrams = new Runnable() {
    			public void run() {
    				for (int j=k;j<wordCounts.size();j+=threadCount) {
    					HashSet<String> batchMatches = new HashSet<String>();
        				tests.add(CheckPermutation(wordCounts.get(j),wordLengths,batchMatches));
        				matches.addAll(batchMatches);
    				}
    			}
    		};
    		var worker = new Thread(workers,GetSinglePermutatonAnagrams);
    		
    		//run thread
    		worker.start();
    	}
		
		//start the status writer
		new Thread(statusWriter).start();
        
        //await threads termination
		while (workers.activeCount()>0)
		try {Thread.sleep(40);} 
		catch (InterruptedException e) {e.printStackTrace();}
	}
	
	
	//checks a single variation from coinchange 
	private static long CheckPermutation(int[] WordLengthCounts, int[] WordLengths, HashSet<String> matches) {
		//initialize ALOT of stuff.
		long testsMade = 0;																//tracks the amount of combinations checked
		int wordSizesCount = 0; for (int i:WordLengthCounts) if (i>0) wordSizesCount++;	//word arrays only containing relevant word lengths
		int[] wordLengthCounts = new int[wordSizesCount],
			  wordLengths = new int[wordSizesCount];
		int index = 0; for(int i=0;i<WordLengthCounts.length;i++) {
			if (WordLengthCounts[i]>0) {
				wordLengthCounts[index] = WordLengthCounts[i];
				wordLengths[index] = WordLengths[i];
				index++;
			}
		}
		int[][] indexes = new int[wordSizesCount][];									//tracks indexes in the word arrays
		for(int i=0;i<indexes.length;i++) indexes[i] = new int[wordLengthCounts[i]];
		StringBuilder builder = new StringBuilder();									//to get words (StringBuilder is faster than String +=)
		
		do {
			//get words factorization
			long result = 1;
			for (int i=0;i<wordLengthCounts.length;i++) {
				for(int j=0;j<wordLengthCounts[i];j++) {
					result *= wordFactorizationByLengthArray[wordLengths[i]][indexes[i][j]];
				}
			}
			
			//check if word is an anagram by doing prime factorization
			testsMade++;
			if (result==anagramFactorization) {
				
				//build anagram and add to list
				matches.add(BuildWords(wordLengthCounts,indexes,wordLengths,builder," "));
			}
		} 
		//increment indexes
		while (IncrementWordIndexes(wordLengthCounts,indexes,wordLengths));
		return testsMade;
	}
	
	//increment indexes. To prevent the list being partial to order, indexes for the same-length words are stored in an array, where when one of these indexes overflow from the word length counts, the index is set to the left index (after incrementing)
	//this creates a list of words impartial to order, as no non unique word combinations are made. This saves significant computation, especially with more same-length words.
	//Basically, for two indexes where the word length count is 3:
	//	[0,0] --> [0,1] --> [0,2] -->
	//	[1,1] --> [1,2] -->
	//	[2,2]
	//Notice that all of these combinations are unique.
	private static boolean IncrementWordIndexes(int[] wordLengthCounts, int[][] indexes, int[] wordLengths) {
		boolean allOverflowed;
		for (int i,j=indexes.length-1;j>=0;j--) {
			i = wordLengthCounts[j];
			
			//go backwards and increment all indexes
			while (--i>=0 && ++indexes[j][i]>=wordByLengthArray[wordLengths[j]].length);
			allOverflowed = (i<0); 
			
			//fix overflows
			while (++i<wordLengthCounts[j]) {
				indexes[j][i] = (i>0? indexes[j][i-1]:0);
			}
			
			if (!allOverflowed) {
				return true;
			}
		}
		//if everything overflowed, then we have reached the end of all indexes.
		return false;
	}
	
	//gets sentence by word indexes
	private static String BuildWords(int[] wordLengthCounts, int[][] indexes, int[] wordLengths, StringBuilder builder, String seperator) {
		//create word
		for (int i=0;i<wordLengthCounts.length;i++) {
			for(int j=0;j<wordLengthCounts[i];j++) {
				builder.append(wordByLengthArray[wordLengths[i]][indexes[i][j]]);
				builder.append(" ");
			}
		}
		String result = builder.toString().substring(0,builder.toString().length()-seperator.length());
		builder.setLength(0);
		return result;
	}
	
	//gets all unique permutations from a list of single permutations (eg. all unique anagrams/anagrams partial to word order from a list of anagrams impartial to word order)
	//we do list of matches so hashset actually makes sense (since Θ(1) insert)
	//this is just Heap's algorithm to keep stuff dynamic (which was the initial goal)
	private static void GetPermutationsOfMatch(HashSet<String> matches) {
		int i;
		for (String match: matches) {
			String[] words = match.split(" ");
			
			//get all permutations of combinations of words
			int[] indexes = new int[words.length];
			i=0;
			while (i<words.length) {
			    if (indexes[i] < i) {
			        Swap(words, i % 2 == 0 ?  0: indexes[i], i);
			        //add permutation
			        matches.add(String.join(" ",words));
			        
			        indexes[i]++;
			        i=0;
			    }
			    else {
			        indexes[i] = 0;
			        i++;
			    }
			}
		}
    }
	
	//swaps two elements
	private static void Swap(int[] elements, int index1, int index2) {
		var temp = elements[index1];
		elements[index1] = elements[index2];
		elements[index2] = temp;
	}
	private static void Swap(String[] elements, int index1, int index2) {
		var temp = elements[index1];
		elements[index1] = elements[index2];
		elements[index2] = temp;
	}
}