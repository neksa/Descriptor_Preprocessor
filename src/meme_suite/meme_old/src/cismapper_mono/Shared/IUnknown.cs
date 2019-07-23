//--------------------------------------------------------------------------------
// File: IUnknown.cs
// Author: Timothy O'Connor
// © Copyright University of Queensland, 2012-2014. All rights reserved.
// License: 
//--------------------------------------------------------------------------------

namespace Shared
{
    using System;

    /// <summary>
    /// Core interface for determining subclass type from a base class
    /// </summary>
    public class IUnknown
    {
        /// <summary>
        /// Initializes a new instance of the <see cref="Shared.IUnknown"/> class.
        /// </summary>
		protected IUnknown()
        {
        }

        /// <summary>
        /// Queries the interface.
        /// </summary>
        /// <returns>The interface.</returns>
        /// <param name="unk">Unk.</param>
        /// <typeparam name="T">The 1st type parameter.</typeparam>
		public static T QueryInterface<T>(IUnknown unk)
		{
			return (T)unk.QueryInterface(typeof(T));
		}

        /// <summary>
        /// Queries the interface.
        /// </summary>
        /// <returns>The interface.</returns>
        /// <param name="t">T.</param>
		virtual protected object QueryInterface(Type t)
		{
			return null;
		}
    }
}

