/**
 * Project: GNormPlus
 * Function: Dictionary lookup by Prefix Tree
 */

package GNormPluslib;

import java.io.*;
import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class PrefixTree
{
	private Tree Tr=new Tree();
	
	/*
	 * Read Dictionary and insert Mention into the Prefix Tree
	 */
	public static HashMap<String, String> StopWord_hash = new HashMap<String, String>();
	
	public void Hash2Tree(HashMap<String, String> ID2Names)
	{
		for(String ID : ID2Names.keySet())  
		{
			String NameColumn[]=ID2Names.get(ID).split("\\|");
			for(int i=0;i<NameColumn.length;i++)
			{
				Tr.insertMention(NameColumn[i],ID);
			}
		}
	}
	public void Dictionary2Tree_Combine(String Filename,String StopWords,String MentionType)	
	{
		try 
		{
			//System.out.println("Dictionary2Tree_Combine : " + Filename);
			
			/** Stop Word */
			BufferedReader br = new BufferedReader(new InputStreamReader(new FileInputStream(StopWords), "UTF-8"));
			String line="";
			while ((line = br.readLine()) != null)  
			{
				StopWord_hash.put(line, "StopWord");
			}
			br.close();	
			
			BufferedReader inputfile = new BufferedReader(new InputStreamReader(new FileInputStream(Filename), "UTF-8"));
			line="";
			//int count=0;
			while ((line = inputfile.readLine()) != null)  
			{
				//count++;
				//if(count%10000==0){	System.out.println(count);	}
				String Column[]=line.split("\t");
				if(Column.length>1)
				{
					Column[0]=Column[0].replace("species:ncbi:","");
					Column[1]=Column[1].replaceAll(" strain=", " ");
					Column[1]=Column[1].replaceAll("[\\W\\-\\_](str\\.|strain|substr\\.|substrain|var\\.|variant|subsp\\.|subspecies|pv\\.|pathovars|pathovar|br\\.|biovar)[\\W\\-\\_]", " ");
					Column[1]=Column[1].replaceAll("[\\(\\)]", " ");
					String SpNameColumn[]=Column[1].split("\\|");
					for(int i=0;i<SpNameColumn.length;i++)
					{
						String tmp = SpNameColumn[i];
						tmp=tmp.replaceAll("[\\W\\-\\_]", "");
						
						/*
						 * Criteria for Species
						 */
						if(	MentionType.equals("Species") &&
							(!SpNameColumn[i].substring(0, 1).matches("[\\W\\-\\_]")) &&
							(!SpNameColumn[i].matches("a[\\W\\-\\_].*")) &&
							tmp.length()>=3
							)
						{
							boolean stopword_boolean=false;
							for(String stopword_RegEx : StopWord_hash.keySet())
							{
								Pattern ptmp = Pattern.compile("^"+stopword_RegEx+"$");
								Matcher mtmp = ptmp.matcher(SpNameColumn[i].toLowerCase());
								if(mtmp.find())
								{
									stopword_boolean=true;
								}
							}
							if(stopword_boolean == false)
							{
								Tr.insertMention(SpNameColumn[i],Column[0]);
							}
						}
						/*
						 * Criteria for Gene
						 */
						else if (MentionType.equals("Gene") &&  
								(!SpNameColumn[i].substring(0, 1).matches("[\\W\\-\\_]")) &&
								tmp.length()>=3
								)
						{
							if(!StopWord_hash.containsKey(SpNameColumn[i].toLowerCase()))
							{
								Tr.insertMention(SpNameColumn[i],Column[0]);
							}
						}
						/*
						 * Criteria for Cell
						 */
						else if (MentionType.equals("Cell") && 
								(!SpNameColumn[i].substring(0, 1).matches("[\\W\\-\\_]")) &&
								tmp.length()>=3
								)
						{
							if(!StopWord_hash.containsKey(SpNameColumn[i].toLowerCase()))
							{
								Tr.insertMention(SpNameColumn[i],Column[0]);
							}
						}
						/*
						 * others
						 */
						else if ((!SpNameColumn[i].substring(0, 1).matches("[\\W\\-\\_]")) &&
								tmp.length()>=3
								)
						{
							if(!StopWord_hash.containsKey(SpNameColumn[i].toLowerCase()))
							{
								Tr.insertMention(SpNameColumn[i],Column[0]);
							}
						}
					}
				}
			}
			inputfile.close();	
		}
		catch(IOException e1){ System.out.println("[Dictionary2Tree_Combine]: Input file is not exist.");}
	}
	public void Dictionary2Tree_UniqueGene(String Filename,String StopWords,String Preifx)	
	{
		try 
		{
			//System.out.println("Dictionary2Tree_UniqueGene : " + Filename);
			
			/** Stop Word */
			BufferedReader br = new BufferedReader(new InputStreamReader(new FileInputStream(StopWords), "UTF-8"));
			String line="";
			while ((line = br.readLine()) != null)  
			{
				StopWord_hash.put(line, "StopWord");
			}
			br.close();	
			
			BufferedReader inputfile = new BufferedReader(new InputStreamReader(new FileInputStream(Filename), "UTF-8"));
			line="";
			//int count=0;
			while ((line = inputfile.readLine()) != null)  
			{
				//count++;
				//if(count%10000==0){	System.out.println(count);	}
				String Column[]=line.split("\t");
				if(Column.length>1)
				{
					if(!StopWord_hash.containsKey(Column[0].toLowerCase()))
					{
						if(Preifx.equals(""))
						{
							Tr.insertMention(Column[0],Column[1]);
						}
						else if(Preifx.equals("Num") && Column[0].matches("[0-9].*"))
						{
							Tr.insertMention(Column[0],Column[1]);
						}
						else if(Preifx.equals("AZNum") && Column[0].matches("[a-z][0-9].*"))
						{
							Tr.insertMention(Column[0],Column[1]);
						}
						else if(Preifx.equals("lo") && Column[0].length()>2 && Column[0].substring(0,2).equals(Preifx))
						{
							if( ! Column[0].matches("loc[0-9]+"))
							{
								Tr.insertMention(Column[0],Column[1]);
							}
						}
						else if(Preifx.equals("un") && Column[0].length()>2 && Column[0].substring(0,2).equals(Preifx))
						{
							if(Column[0].length()>=6 && Column[0].substring(0,6).equals("unchar"))
							{
								// remove uncharacterized
							}
							else
							{
								Tr.insertMention(Column[0],Column[1]);
							}
						}
						else if(Column[0].length()>2 && Column[0].substring(0,2).equals(Preifx))
						{
							Tr.insertMention(Column[0],Column[1]);
						}
					}
				}
			}
			inputfile.close();	
		}
		catch(IOException e1){ System.out.println("[Dictionary2Tree_UniqueGene]: Input file is not exist.");}
	}
	public void Dictionary2Tree_UniqueSpecies(String Filename,String StopWords,String Preifx)	
	{
		try 
		{
			//System.out.println("Dictionary2Tree_UniqueGene : " + Filename);
			
			/** Stop Word */
			BufferedReader br = new BufferedReader(new InputStreamReader(new FileInputStream(StopWords), "UTF-8"));
			String line="";
			while ((line = br.readLine()) != null)  
			{
				StopWord_hash.put(line, "StopWord");
			}
			br.close();	
			
			BufferedReader inputfile = new BufferedReader(new InputStreamReader(new FileInputStream(Filename), "UTF-8"));
			line="";
			while ((line = inputfile.readLine()) != null)  
			{
				//count++;
				//if(count%10000==0){	System.out.println(count);	}
				String Column[]=line.split("\t");
				if(Column.length>1)
				{
					if(!StopWord_hash.containsKey(Column[0].toLowerCase()))
					{
						if(Preifx.equals("")) //all
						{
							if(Column[0].matches(".*[\\W\\-\\_](str\\.|strain|substr\\.|substrain|var\\.|variety|variant|subsp\\.|subspecies|pv\\.|pathovars|pathovar|br\\.|biovar)[\\W\\-\\_].*"))
							{
								String mention_rev=Column[0].replaceAll("[\\W\\-\\_](str\\.|strain|substr\\.|substrain|var\\.|variety|variant|subsp\\.|subspecies|pv\\.|pathovars|pathovar|br\\.|biovar)[\\W\\-\\_]", " ");
								String mention_tmp=mention_rev.replaceAll("[\\W\\-\\_]","");
								if(mention_tmp.length()>=10)
								{
									Tr.insertMention(mention_rev,Column[1]);
								}
							}
							else
							{
								Tr.insertMention(Column[0],Column[1]); // mention, id
							}
							
						}
						else if(Column[0].matches("[0-9][0-9].*"))
						{
							if(Preifx.equals("Num"))
							{
								if(Column[0].matches(".*[\\W\\-\\_](str\\.|strain|substr\\.|substrain|var\\.|variety|variant|subsp\\.|subspecies|pv\\.|pathovars|pathovar|br\\.|biovar)[\\W\\-\\_].*"))
								{
									String mention_rev=Column[0].replaceAll("[\\W\\-\\_](str\\.|strain|substr\\.|substrain|var\\.|variety|variant|subsp\\.|subspecies|pv\\.|pathovars|pathovar|br\\.|biovar)[\\W\\-\\_]", " ");
									String mention_tmp=mention_rev.replaceAll("[\\W\\-\\_]","");
									if(mention_tmp.length()>=10)
									{
										Tr.insertMention(mention_rev,Column[1]);
									}
								}
								else
								{
									Tr.insertMention(Column[0],Column[1]); // mention, id
								}
							}
						}
						/*
						else if(Column[0].matches("[a-z][0-9].*"))
						{
							if(Preifx.equals("AZNum"))
							{
								if(Column[0].matches(".*[\\W\\-\\_](str\\.|strain|substr\\.|substrain|var\\.|variety|variant|subsp\\.|subspecies|pv\\.|pathovars|pathovar|br\\.|biovar)[\\W\\-\\_].*"))
								{
									String mention_rev=Column[0].replaceAll("[\\W\\-\\_](str\\.|strain|substr\\.|substrain|var\\.|variety|variant|subsp\\.|subspecies|pv\\.|pathovars|pathovar|br\\.|biovar)[\\W\\-\\_]", " ");
									String mention_tmp=mention_rev.replaceAll("[\\W\\-\\_]","");
									if(mention_tmp.length()>=10)
									{
										Tr.insertMention(mention_rev,Column[1]);
									}
								}
								else
								{
									Tr.insertMention(Column[0],Column[1]); // mention, id
								}
							}
						}
						*/
						else if(Column[0].matches("[a-z][a-z].*"))
						{
							if(Column[0].length()>2 && Column[0].substring(0,2).equals(Preifx))
							{
								if(Column[0].matches(".*[\\W\\-\\_](str\\.|strain|substr\\.|substrain|var\\.|variety|variant|subsp\\.|subspecies|pv\\.|pathovars|pathovar|br\\.|biovar)[\\W\\-\\_].*"))
								{
									String mention_rev=Column[0].replaceAll("[\\W\\-\\_](str\\.|strain|substr\\.|substrain|var\\.|variety|variant|subsp\\.|subspecies|pv\\.|pathovars|pathovar|br\\.|biovar)[\\W\\-\\_]", " ");
									String mention_tmp=mention_rev.replaceAll("[\\W\\-\\_]","");
									if(mention_tmp.length()>=10)
									{
										Tr.insertMention(mention_rev,Column[1]);
									}
								}
								else
								{
									Tr.insertMention(Column[0],Column[1]); // mention, id
								}
							}
						}
						else if(Preifx.equals("Others"))
						{
							if(Column[0].matches(".*[\\W\\-\\_](str\\.|strain|substr\\.|substrain|var\\.|variety|variant|subsp\\.|subspecies|pv\\.|pathovars|pathovar|br\\.|biovar)[\\W\\-\\_].*"))
							{
								String mention_rev=Column[0].replaceAll("[\\W\\-\\_](str\\.|strain|substr\\.|substrain|var\\.|variety|variant|subsp\\.|subspecies|pv\\.|pathovars|pathovar|br\\.|biovar)[\\W\\-\\_]", " ");
								String mention_tmp=mention_rev.replaceAll("[\\W\\-\\_]","");
								if(mention_tmp.length()>=10)
								{
									Tr.insertMention(mention_rev,Column[1]);
								}
							}
							else
							{
								Tr.insertMention(Column[0],Column[1]); // mention, id
							}
						}
					}
				}
			}
			inputfile.close();	
		}
		catch(IOException e1){ System.out.println("[Dictionary2Tree_UniqueGene]: Input file is not exist.");}
	}
	public void TreeFile2Tree(String Filename)	
	{
		try 
		{
			//System.out.println("TreeFile2Tree : " + Filename);
			
			BufferedReader inputfile = new BufferedReader(new InputStreamReader(new FileInputStream(Filename), "UTF-8"));
			String line="";
			int count=0;
			while ((line = inputfile.readLine()) != null)  
			{
				String Anno[]=line.split("\t");
				if(Anno.length<2){System.out.println(count+"\t"+line);} //check error
				String LocationInTree = Anno[0];
				String token = Anno[1];
				String identifier="";
				if(Anno.length==3)
				{
					identifier = Anno[2];
				}
				String LocationsInTree[]=LocationInTree.split("-");
				TreeNode tmp = Tr.root;
				for(int i=0;i<LocationsInTree.length-1;i++)
				{
					tmp=tmp.links.get(Integer.parseInt(LocationsInTree[i])-1);
				}
				tmp.InsertToken(token,identifier);
				//if(count%10000==0){System.out.println(count);}
				count++;
			}
			inputfile.close();	
		}
		catch(IOException e1){ System.out.println("[TreeFile2Tee]: Input file: "+ Filename +" is not exist.");}
	}
	
	/*
	 * Search target mention in the Prefix Tree
	 */
	public String MentionMatch(String Mentions)
	{
		ArrayList<String> location = new ArrayList<String>();
		String Menlist[]=Mentions.split("\\|");
		for(int m=0;m<Menlist.length;m++)
		{
			String Mention=Menlist[m];
			String Mention_lc=Mention.toLowerCase();
			Mention_lc = Mention_lc.replaceAll("[\\W\\-\\_]+", "");
			Mention_lc = Mention_lc.replaceAll("([0-9])([a-z])", "$1 $2");
			Mention_lc = Mention_lc.replaceAll("([a-z])([0-9])", "$1 $2");
			String Tkns[]=Mention_lc.split(" ");
			
			int PrefixTranslation=0;
			int i=0;
			boolean find=false;
			TreeNode tmp = Tr.root;
			
			while( i<Tkns.length && tmp.CheckChild(Tkns[i],PrefixTranslation)>=0) //Find Tokens in the links
			{
				if(i == Tkns.length-1){PrefixTranslation = 1;}
				tmp=tmp.links.get(tmp.CheckChild(Tkns[i],PrefixTranslation)); //move point to the link
				find=true;
				i++;
			}
			if(find == true)
			{
				if(i==Tkns.length)
				{
					if(!tmp.Concept.equals(""))
					{
						return tmp.Concept;
					}
					else
					{
						return "-1";
						//gene id is not found.
					}
				}
				else
				{
					return "-2";
					//the gene mention matched a substring in PrefixTree.
				}
			}
			else
			{
				return "-3";
				//mention is not found
			}
		}
		return "-3"; //mention is not found
	}
	
	/*
	 * Search target mention in the Prefix Tree
	 */
	public String MentionMatch_species(String Mentions)
	{
		ArrayList<String> location = new ArrayList<String>();
		String Menlist[]=Mentions.split("\\|");
		for(int m=0;m<Menlist.length;m++)
		{
			String Mention=Menlist[m];
			String Mention_lc=Mention.toLowerCase();
			Mention_lc = Mention_lc.replaceAll("[\\W\\-\\_]+", " ");
			Mention_lc = Mention_lc.replaceAll("([0-9])([a-z])", "$1 $2");
			Mention_lc = Mention_lc.replaceAll("([a-z])([0-9])", "$1 $2");
			Mention_lc = Mention_lc.replaceAll("^[ ]+", "");
			Mention_lc = Mention_lc.replaceAll("[ ]+$", "");
			String Tkns[]=Mention_lc.split(" ");
			
			int PrefixTranslation=0;
			int i=0;
			boolean find=false;
			TreeNode tmp = Tr.root;
			
			while( i<Tkns.length && tmp.CheckChild(Tkns[i],PrefixTranslation)>=0) //Find Tokens in the links
			{
				if(i == Tkns.length-1){PrefixTranslation = 1;}
				tmp=tmp.links.get(tmp.CheckChild(Tkns[i],PrefixTranslation)); //move point to the link
				find=true;
				i++;
			}
			if(find == true)
			{
				if(i==Tkns.length)
				{
					if(!tmp.Concept.equals(""))
					{
						return tmp.Concept;
					}
					else
					{
						return "-1";
						//gene id is not found.
					}
				}
				else
				{
					return "-2";
					//the gene mention matched a substring in PrefixTree.
				}
			}
			else
			{
				return "-3";
				//mention is not found
			}
		}
		return "-3"; //mention is not found
	}
	
	/*
	 * Search target mention in the Prefix Tree
	 * ConceptType: Species|Genus|Cell|CTDGene
	 */
	public ArrayList<String> SearchMentionLocation(String Doc,String ConceptType)
	{
		ArrayList<String> location = new ArrayList<String>();
		Doc=Doc+" XXXX XXXX";
		String Doc_org=Doc;
		Doc=Doc.toLowerCase();
		String Doc_lc=Doc;
		Doc = Doc.replaceAll("([0-9])([A-Za-z])", "$1 $2");
		Doc = Doc.replaceAll("([A-Za-z])([0-9])", "$1 $2");
		Doc = Doc.replaceAll("[\\W^;:,]+", " ");
		
		/* = keep special characters =
		 * 
		String regex="\\s+|(?=\\p{Punct})|(?<=\\p{Punct})";
		String DocTkns[]=Doc.split(regex);
		 */
		
		String DocTkns[]=Doc.split(" ");
		int Offset=0;
		int Start=0;
		int Last=0;
		int FirstTime=0;
		
		while(Doc_lc.length()>0 && Doc_lc.substring(0,1).matches("[\\W]")) //clean the forward whitespace
		{
			Doc_lc=Doc_lc.substring(1);
			Offset++;
		}
		
		for(int i=0;i<DocTkns.length;i++)
		{
			//System.out.println(i+"\t"+Start+"\t"+Last+"\t"+Offset+"\t"+Doc_lc);
			
			int pre_i=i;
			int pre_Start=Start;
			int pre_Last=Last;
			String pre_Doc_lc=Doc_lc;
			int pre_Offset=Offset;
			
			TreeNode tmp = Tr.root;
			boolean find=false;
			int PrefixTranslation=2;
			if(ConceptType.equals("Species"))
			{
				PrefixTranslation=3;
			}
			int ConceptFound=i; //Keep found concept
			String ConceptFound_STR="";//Keep found concept
			int FirstTime_while = -1;
			
			while( tmp.CheckChild(DocTkns[i],PrefixTranslation)>=0 ) //Find Tokens in the links
			{
				FirstTime_while++;
				tmp=tmp.links.get(tmp.CheckChild(DocTkns[i],PrefixTranslation)); //move point to the link
				if(Start==0 && FirstTime>0){Start = Offset;} //Start <- Offset 
				if(Doc_lc.length()>=DocTkns[i].length() && Doc_lc.substring(0,DocTkns[i].length()).equals(DocTkns[i]))
				{
					if(DocTkns[i].length()>0)
					{
						Doc_lc=Doc_lc.substring(DocTkns[i].length());
						Offset=Offset+DocTkns[i].length();
					}
				}
				Last = Offset;
				while(Doc_lc.length()>0 && Doc_lc.substring(0,1).matches("[\\W]")) //clean the forward whitespace
				{
					Doc_lc=Doc_lc.substring(1);
					Offset++;
				}
				i++;
				
				if(ConceptType.equals("Species"))
				{
					if(i<DocTkns.length-3 && DocTkns[i].matches("(str|strain|substr|substrain|subspecies|subsp|var|variant|pathovars|pv|biovar|bv)"))
					{
						Doc_lc=Doc_lc.substring(DocTkns[i].length());
						Offset=Offset+DocTkns[i].length();
						Last = Offset;
						while(Doc_lc.length()>0 && Doc_lc.substring(0,1).matches("[\\W]")) //clean the forward whitespace
						{
							Doc_lc=Doc_lc.substring(1);
							Offset++;
						}
						i++;
					}
				}
				
				if(!tmp.Concept.equals("") && (Last-Start>0)) //Keep found concept
				{
					if(Last<Doc_org.length())
					{
						ConceptFound=i;
						ConceptFound_STR=Start+"\t"+Last+"\t"+Doc_org.substring(Start, Last)+"\t"+tmp.Concept;
						//System.out.println(ConceptFound_STR);
					}
				}
				
				find=true;
				if(i>=DocTkns.length){break;}
				else if(i==DocTkns.length-1){PrefixTranslation=2;}
				
				//System.out.println(i+"\t"+Start+"\t"+Last+"\t("+FirstTime_while+")\t"+Offset+"\t"+Doc_lc);
				
				if(FirstTime_while==0) // first matched token
				{
					pre_i=i;
					pre_Start=Start;
					pre_Last=Last;
					pre_Doc_lc=Doc_lc;
					pre_Offset=Offset;
				}
			}
			
			if(find == true)
			{
				//System.out.println(find+"\t"+FirstTime_while+"\t"+Start+"\t"+Last+"\t"+Doc_org.substring(Start, Last)+"\t"+tmp.Concept);
				if(!tmp.Concept.equals("")) //the last matched token has concept id 
				{
					if(Last<Doc_org.length() && Last>Start)
					{
						location.add(Start+"\t"+Last+"\t"+Doc_org.substring(Start, Last)+"\t"+tmp.Concept);
					}
				}
				else
				{
					if(!ConceptFound_STR.equals("")) //Keep found concept
					{
						location.add(ConceptFound_STR);
						i = ConceptFound + 1;
					}
					
					if(FirstTime_while>=1)
					{
						i=pre_i;
						Start=pre_Start;
						Last=pre_Last;
						Doc_lc=pre_Doc_lc;
						Offset=pre_Offset;
					}
				}
				Start=0;
				Last=0;
				if(i>0){i--;}
				ConceptFound=i; //Keep found concept
				ConceptFound_STR="";//Keep found concept
			}
			else //if(find == false)
			{
				//System.out.println(find+"\t"+FirstTime_while+"\t"+Start+"\t"+Last+"\t"+Doc_org.substring(Start, Last)+"\t"+tmp.Concept);
				
				if(FirstTime_while>=1 && tmp.Concept.equals(""))
				{
					i=pre_i;
					Start=pre_Start;
					Last=pre_Last;
					Doc_lc=pre_Doc_lc;
					Offset=pre_Offset;
				}
				
				if(Doc_lc.length()>=DocTkns[i].length() && Doc_lc.substring(0,DocTkns[i].length()).equals(DocTkns[i]))
				{
					if(DocTkns[i].length()>0)
					{
						Doc_lc=Doc_lc.substring(DocTkns[i].length());
						Offset=Offset+DocTkns[i].length();
					}
				}
			}
			
			while(Doc_lc.length()>0 && Doc_lc.substring(0,1).matches("[\\W]")) //clean the forward whitespace
			{
				Doc_lc=Doc_lc.substring(1);
				Offset++;
			}
			FirstTime++;
			
			//System.out.println();
		}
		return location;
	}
	
	/*
	 * Print out the Prefix Tree
	 */
	public String PrintTree()
	{
		return Tr.PrintTree_preorder(Tr.root,"");
	}
	
	public void SaveTree(String outputfile) throws IOException
	{
		BufferedWriter fr = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(outputfile), "UTF-8"));
		Tr.SaveTree_preorder(Tr.root,"",fr);
		fr.close();
	}
	
	
	public void insertMention(String Mention, String Identifier)
	{
		Tr.insertMention(Mention,Identifier);
	}
}

class Tree 
{
	/*
	 * Prefix Tree - root node
	 */
	public TreeNode root;
	
	public Tree() 
	{ 
		root = new TreeNode("-ROOT-"); 
	}
	
	/*
	 * Insert mention into the tree
	 */
	public void insertMention(String Mention, String Identifier)
	{
		Mention=Mention.toLowerCase();
		
		Mention = Mention.replaceAll("([0-9])([A-Za-z])", "$1 $2");
		Mention = Mention.replaceAll("([A-Za-z])([0-9])", "$1 $2");
		Mention = Mention.replaceAll("[\\W\\-\\_]+", " ");
		/* = keep special characters =
		 * 
		String regex="\\s+|(?=\\p{Punct})|(?<=\\p{Punct})";
		String Tokens[]=Mention.split(regex);
		 */
		String Tokens[]=Mention.split(" ");
		TreeNode tmp = root;
		for(int i=0;i<Tokens.length;i++)
		{
			if(tmp.CheckChild(Tokens[i],0)>=0)
			{
				tmp=tmp.links.get( tmp.CheckChild(Tokens[i],0) ); //go through next generation (exist node)
				if(i == Tokens.length-1)
				{
					tmp.Concept=Identifier;
				}
			}
			else //not exist
			{
				if(i == Tokens.length-1)
				{
					tmp.InsertToken(Tokens[i],Identifier);
				}
				else
				{
					tmp.InsertToken(Tokens[i]);
				}
				tmp=tmp.links.get(tmp.NumOflinks-1); //go to the next generation (new node)
			}
		}
	}
	
	/*
	 * Print the tree by pre-order
	 */
	public String PrintTree_preorder(TreeNode node, String LocationInTree)
	{
		String opt="";
		if(!node.token.equals("-ROOT-"))//Ignore root
		{
			if(node.Concept.equals(""))
			{
				opt=opt+LocationInTree+"\t"+node.token+"\n";
			}
			else
			{
				opt=opt+LocationInTree+"\t"+node.token+"\t"+node.Concept+"\n";
			}
		} 
		if(!LocationInTree.equals("")){LocationInTree=LocationInTree+"-";}
		for(int i=0;i<node.NumOflinks;i++)
		{
			opt=opt+PrintTree_preorder(node.links.get(i),LocationInTree+(i+1));
		}
		return opt;
	}
	
	/*
	 * Print the tree by pre-order
	 */
	public void SaveTree_preorder(TreeNode node, String LocationInTree, BufferedWriter fr) throws IOException
	{
		if(!node.token.equals("-ROOT-"))//Ignore root
		{
			if(node.Concept.equals(""))
			{
				fr.write(LocationInTree+"\t"+node.token+"\n");
			}
			else
			{
				fr.write(LocationInTree+"\t"+node.token+"\t"+node.Concept+"\n");
			}
		} 
		if(!LocationInTree.equals("")){LocationInTree=LocationInTree+"-";}
		for(int i=0;i<node.NumOflinks;i++)
		{
			SaveTree_preorder(node.links.get(i),LocationInTree+(i+1),fr);
		}
	}
}

class TreeNode 
{
	String token; //token of the node
	int NumOflinks; //Number of links
	public String Concept;
	HashMap<String,Integer> Hashs;
	ArrayList<TreeNode> links;
	
	public TreeNode(String Tok,String ID)
	{
		token = Tok;
		NumOflinks = 0;
		Concept = ID;
		links = new ArrayList<TreeNode>();/*link*/
		Hashs = new HashMap<String,Integer>();/*hash*/
	}
	public TreeNode(String Tok)
	{
		token = Tok;
		NumOflinks = 0;
		Concept = "";
		links = new ArrayList<TreeNode>();/*link*/
		Hashs = new HashMap<String,Integer>();/*hash*/
	}
	public TreeNode()
	{
		token = "";
		NumOflinks = 0;
		Concept = "";
		links = new ArrayList<TreeNode>();/*link*/
		Hashs = new HashMap<String,Integer>();/*hash*/
	}
	
	public String toString()
	{
		return (token+"\t"+Concept);
	}
	
	/*
	 * Insert an new node under the target node
	 */
	public void InsertToken(String Tok)
	{
		TreeNode NewNode = new TreeNode(Tok);
		
		/*link*/
		links.add(NewNode);
		
		/*hash*/
		Hashs.put(Tok, NumOflinks);
		
		NumOflinks++;
	}
	public void InsertToken(String Tok,String ID)
	{
		TreeNode NewNode = new TreeNode(Tok,ID);
		/*link*/
		links.add(NewNode);
		
		/*hash*/
		Hashs.put(Tok, NumOflinks);
		
		NumOflinks++;
	}
	
	/*
	 * Check the tokens of children
	 */
	public int CheckChild(String Tok, Integer PrefixTranslation)
	{
		if(Hashs.containsKey(Tok))
		{
			return(Hashs.get(Tok));
		}
		
		if(PrefixTranslation == 1 && Tok.matches("(alpha|beta|gamam|[abg]|[12])")) // SuffixTranslationMap
		{
			if(Hashs.containsKey(GNormPlus.SuffixTranslationMap_hash.get(Tok)))
			{
				return(Hashs.get(GNormPlus.SuffixTranslationMap_hash.get(Tok)));
			}
			
		}
		else if(PrefixTranslation == 2 && Tok.matches("[1-5]")) // for CTDGene feature
		{
			for(int i=0;i<links.size();i++)
	        {
				if(links.get(i).token.matches("[1-5]"))
				{
					return(i);
				}
	        }
			
			for(int i=1;i<=5;i++)
			{
				if(Hashs.containsKey(i)){return(Hashs.get(i));}
			}
		}
		
		return(-1);
	}
}
	