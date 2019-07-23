
namespace Genomics
{
    using System;
    using System.Collections.Generic;
    using System.Linq;

    public class TrieNode<AlphabetType> where AlphabetType : Alphabet
    {
        public TrieNode()
        {
            this.Children = new Dictionary<char, TrieNode<AlphabetType>>();
        }

        public char Letter { get; set; }

        public void Add(string s)
        {
            if (s.Length == 0)
            {
                this.EndCount++;
                return;
            }

            var prefix = s.First();
            var suffix = s.Substring(1, s.Length - 1);

            if (!this.Children.ContainsKey(prefix))
            {
                this.Children.Add(prefix, new TrieNode<AlphabetType>());
            }

            this.Children[prefix].Add(suffix);
        }

        public int EndCount { get; set; }

        public Dictionary<char, TrieNode<AlphabetType>> Children { get; set; }
    }
}

