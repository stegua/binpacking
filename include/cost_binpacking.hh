#ifndef __MY_PACK_PROPGATOR
#define __MY_PACK_PROPGATOR

#include <gecode/int.hh>

/**
 * \namespace Gecode::Int::CostBinPacking
 * \brief %Bin-packing propagators
 */

namespace Gecode { namespace Int { namespace CostBinPacking {

  /**
   * \brief Item combining bin and size information
   */
  class Item : public DerivedView<IntView> {
  protected:
    using DerivedView<IntView>::x;
    /// Size of item
    int s;
  public:
    /// Default constructor
    Item(void);
    /// Constructor
    Item(IntView b, int s);

    /// Return bin of item
    IntView bin(void) const;
    /// Set bin of item to \a b
    void bin(IntView b);
    /// Return size of item
    int size(void) const;
    /// Set size of item to \a s
    void size(int s);

    /// Update item during cloning
    void update(Space& home, bool share, Item& i);
  };

  /// Whether two items are the same
  bool same(const Item& i, const Item& j);
  /// Test whether one item is before another
  bool before(const Item& i, const Item& j);

  /// For sorting according to size
  bool operator <(const Item& i, const Item& j);


  /// Size sets
  class SizeSet {
  protected:
    /// Number of size entries in the set
    int n;
    /// Total size of the set
    int t;
    /// Array of sizes (will have more elements)
    int* s;
  public:
    /// Default constructor
    SizeSet(void);
    /// Initialize for at most \a n_max items
    SizeSet(Region& region, int n_max);
    /// Add new size \a s
    void add(int s);
    /// Return cardinality of set (number of entries)
    int card(void) const;
    /// Return total size
    int total(void) const;
    /// Return size of item \a i
    int operator [](int i) const;
  };

  /// Size sets with one element discarded
  class SizeSetMinusOne : public SizeSet {
  protected:
    /// Position of discarded item
    int p;
  public:
    /// Default constructor
    SizeSetMinusOne(void);
    /// Initialize for at most \n n_max entries
    SizeSetMinusOne(Region& region, int n);
    /// Discard size \a s
    void minus(int s);
    /// Return cardinality of set (number of entries)
    int card(void) const;
    /// Return total size
    int total(void) const;
    /// Return size of item \a i
    int operator [](int i) const;
  };


  /**
   * \brief Bin-packing propagator
   *
   * The algorithm is taken from:
   *   Paul Shaw. A Constraint for Bin Packing. CP 2004.
   *
   * Requires \code #include <gecode/int/bin-packing.hh> \endcode
   *
   * \ingroup FuncIntProp
   */
  class Pack : public Propagator {
  protected:
    /// Views for load of bins
    ViewArray<OffsetView> l;
    /// Items with bin and size
    ViewArray<Item> x;
    /// Cost variable
    IntView z;
    /// Cost of assign processi i to machine j
    IntSharedArray  D;
    /// Constructor for posting
    Pack(Home home, ViewArray<OffsetView>& l, ViewArray<Item>& x, IntView& z, const IntSharedArray& D);
    /// Constructor for cloning \a p
    Pack(Space& home, bool share, Pack& p);
  public:
    /// Post propagator for loads \a l and items \a x
    GECODE_INT_EXPORT 
    static ExecStatus post(Home home, 
                           ViewArray<OffsetView>& l, 
                           ViewArray<Item>&       x, 
                           IntView&               z,
                           const IntSharedArray&  D);
    /// Perform propagation
    GECODE_INT_EXPORT 
    virtual ExecStatus propagate(Space& home, const ModEventDelta& med);
    /// Cost function
    GECODE_INT_EXPORT 
    virtual PropCost cost(const Space& home, const ModEventDelta& med) const;
    /// Copy propagator during cloning
    GECODE_INT_EXPORT 
    virtual Actor* copy(Space& home, bool share);
    /// Destructor
    virtual size_t dispose(Space& home);
  };

}}

GECODE_INT_EXPORT void
cost_binpacking(Home home, 
      const IntVarArgs& l, 
      const IntVarArgs& x, 
      const IntVar&     z,
      const IntArgs&    w, 
      const IntSharedArray& D,
      IntConLevel icl=ICL_DEF);

}

#include "propagate.hh"

#endif


