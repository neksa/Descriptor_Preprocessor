
namespace Genomics
{
    using System;
    using System.Collections.Generic;
    using System.Linq;
    using Shared;
   

    public class Trie<AlphabetType> where AlphabetType : NucleotideAlphabet
    {
        private NucleotideAlphabet alphabet;

        public Trie()
        {
        }

        public void AddSequence(string sequence, int length, int id)
        {
            var reverseSequence = this.Alphabet.ReverseComplement(sequence);

            for (int i = 0; i < sequence.Length - length; i++)
            {
                this.root.Add(sequence.Substring(i, length));
                this.root.Add(reverseSequence.Substring(i, length));
            }
        }

        public NucleotideAlphabet Alphabet
        {
            get
            {
                return Helpers.CheckInit(
                    ref this.alphabet,
                    () => NucleotideAlphabets.Get(typeof(AlphabetType)));
            }
        }

        public int Count()
        {
            return this.Count(root);
        }

        public int NodeCount()
        {
            return this.NodeCount(root);
        }

        public int CharCount()
        {
            return this.CharCount(root, 0);
        }

        private int Count(TrieNode<AlphabetType> node)
        {
            return node.EndCount + node.Children.Sum(x => this.Count(x.Value));
        }

        private int NodeCount(TrieNode<AlphabetType> node)
        {
            return 1 + node.Children.Sum(x => this.NodeCount(x.Value));
        }

        private int CharCount(TrieNode<AlphabetType> node, int depth)
        {
            return depth * node.EndCount + node.Children.Sum(x => this.CharCount(x.Value, depth + 1));
        }

        private TrieNode<AlphabetType> root = new TrieNode<AlphabetType>();
    }
}

