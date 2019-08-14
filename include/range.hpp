/* range.hpp
 *                                                                               
 * Copyright (c) 2013 Alexander Duchene <aduche4@tigers.lsu.edu>                                                     
 *                                                                              
 * This piece of software was created as part of the Drosophila  Population 
 * Genomics Project opensource agreement.
 * 
 * Permission is hereby granted, free of charge, to any person obtaining a copy 
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights 
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell    
 * copies of the Software, and to permit persons to whom the Software is        
 * furnished to do so, subject to the following conditions:                     
 * 
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

#ifndef RANGE_ITERATOR_HPP
#define RANGE_ITERATOR_HPP
#include <iterator>

namespace rangepp {

  template<typename value_t>
  class range_impl{
  private:
    value_t rbegin;
    value_t rend;
    value_t step;
    int step_end;
  public:
    constexpr range_impl(value_t begin, value_t end, value_t step=1):
    rbegin(begin), rend(end), step(step), step_end((rend-rbegin)/step) {
      if(rbegin+step_end*step != rend){
	step_end++;
      }
    }

    class iterator:
      public std::iterator<std::random_access_iterator_tag,value_t>
    {
    private:
      value_t current_value;
      int current_step;
      range_impl& parent;
    public:
      constexpr iterator(int start,range_impl& p): current_value(p.rbegin + start*p.step),
						   current_step(start), parent(p) {
      }
      constexpr value_t operator*() {
	return current_value;
      }
      constexpr const iterator* operator++(){
	current_value+=parent.step;
	current_step++;
	return this;
      }
      constexpr const iterator* operator++(int){
	current_value+=parent.step;
	current_step++;
	return this;
      }
      constexpr bool operator==(const iterator& other) {
	return current_step==other.current_step;
      }
      constexpr bool operator!=(const iterator& other) {
	return current_step!=other.current_step;
      }
      constexpr iterator operator+(int s) {
	iterator ret=*this;
	ret.current_step+=s;
	ret.current_value+=s*parent.step;
	return ret;
      }
      constexpr iterator operator-(int s){
	iterator ret=*this;
	ret.current_step-=s;
	ret.current_value-=s*parent.step;
	return ret;
      }
      constexpr const iterator* operator--(){
	current_value-=parent.step;
	current_step--;
	return this;}
      constexpr iterator operator--(int){
	iterator old=*this;
	current_value-=parent.step;
	current_step--;
	return old;
      }
    };

    constexpr iterator begin(){
      return iterator(0,*this);
    }
    constexpr iterator end(){
      return iterator(step_end,*this);
    }

    constexpr value_t operator[](int s){
      return rbegin+s*step;
    }

    constexpr int size(){
      return step_end;
    }
  };
}
template<typename vt,typename other>
constexpr auto range(other begin, other end, vt stepsize)->rangepp::range_impl<decltype(begin+end+stepsize)>
{
    
  return rangepp::range_impl<decltype(begin+end+stepsize)>(begin,end,stepsize);
}

template<typename b,typename e>
constexpr auto range(b begin, e end) -> rangepp::range_impl<decltype(begin+end)>
{
  return rangepp::range_impl<decltype(begin+end)>(begin,end,1);
}

template<typename e>
constexpr rangepp::range_impl<e> range(e end){
  return rangepp::range_impl<e>(0,end,1);
}

#endif
