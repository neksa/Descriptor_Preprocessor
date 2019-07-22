// <copyright file="AnalysisElement.cs" 
//            company="The University of Queensland"
//            author="Timothy O'Connor">
//     Copyright © The University of Queensland, 2012-2014. All rights reserved.
// </copyright>
// License: 
//--------------------------------------------------------------------------------

namespace Analyses
{
    using System;

    /// <summary>
    /// Analysis element.
    /// </summary>
    abstract public class AnalysisElement<TElement>
    {
        /// <summary>
        /// Gets or sets the property prefix.
        /// </summary>
        /// <value>The property prefix.</value>
        protected string PropertyPrefix { get; set; }

        /// <summary>
        /// The element.
        /// </summary>
        protected TElement element;

        /// <summary>
        /// Initializes a new instance of the <see cref="Analyses.AnalysisElement`1"/> class.
        /// </summary>
        protected AnalysisElement()
        {
            this.PropertyPrefix = string.Empty;
        }

        /// <summary>
        /// Initializes a new instance of the <see cref="Analyses.AnalysisElement`1"/> class.
        /// </summary>
        /// <param name="propertyPrefix">Property prefix.</param>
        protected AnalysisElement(string propertyPrefix)
        {
            this.PropertyPrefix = propertyPrefix;
        }

        /// <summary>
        /// Gets or sets the element.
        /// </summary>
        /// <value>The element.</value>
        public virtual TElement Element
        {
            get
            {
                return this.element;
            }

            set
            {
                this.element = value;
            }
        }

        /// <param name="e">E.</param>
        public static implicit operator TElement(AnalysisElement<TElement> e) 
        {
            return e.Element;
        }

        /// <summary>
        /// Register the specified analysis.
        /// </summary>
        /// <param name="analysis">Analysis.</param>
        public abstract void Register(BaseAnalysis analysis);

        /// <summary>
        /// Registers the property.
        /// </summary>
        /// <param name="analysis">Analysis.</param>
        /// <param name="name">Name.</param>
        /// <param name="action">Action.</param>
        protected void RegisterProperty(BaseAnalysis analysis, string name, Action<string> action)
        {
            analysis.ElementArgRegistry.Add(this.PropertyPrefix + name, action);
        }
    }
}

