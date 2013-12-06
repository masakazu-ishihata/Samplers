#!/usr/bin/env ruby
# -*- coding: utf-8 -*-

########################################
# Array
########################################
class Array
  def sum; self.inject(0){|s, i| s+=i}; end
end

########################################
# Gamma Distribution
########################################
class Gamma
  #### new ####
  attr_reader :k, :th
  def initialize(_k, _th)
    @k  = _k.to_f
    @th = _th.to_f
  end

  #### sample ####
  def sample
    if @k < 1
      begin
        u = rand ** (1 / @k)
        v = rand ** (1 / (1 - @k))
      end while u + v > 1
      x = -Math.log(rand) * (u / (u + v))
    else
      b = @k - 1
      c = 3 * @k - 0.75
      begin
        u = rand
        v = rand
        w = u * (1 - u)
        y = Math.sqrt(c / w) * (u - 0.5)
        x = b + y
        redo if x < 0
        z = 64 * w * w * w * v * v
        break if (z <= 1-(2*y*y)/x) || (Math.log(z) < 2*(b*Math.log(x/b)-y))
      end while true
    end
    @th * x
  end

  #### show ####
  def show; puts "#{@k}, #{@th}"; end
end

########################################
# Gamma Distribution
########################################
class Dirichlet
  #### new ####
  attr_reader :al
  def initialize(_al); @al = _al.clone; end
  def Dirichlet.new_simple(_d, _al)
    Dirichlet.new( Array.new(_d){|d| _al} )
  end

  #### sample ####
  def sample
    th = Array.new(@al.size){|i| Gamma.new(@al[i], 1.0).sample}
    sum = th.sum
    th.map!{|i| i / sum}
  end

  #### show ####
  def show; puts "[#{@al.join(", ")}]"; end
end

########################################
# Categorical 分布
########################################
class Categorical
  #### new ####
  attr_reader :th
  def initialize(_theta)
    @th = _theta.clone
    normalize
  end
  def Categorical.new_simple(_d, _al)
    Categorical.new( Dirichlet.new_simple(_d, _al).sample )
  end

  #### normalize ####
  def normalize
    sum = @th.sum.to_f
    @th.map!{|i| i / sum}
  end

  #### pof ####
  def pof(_i); @th[_i]; end

  #### sample ####
  def sample
    r = rand
    sum = 0
    for i in 0..@th.size-1
      sum += @th[i]
      break if r < sum
    end
    i
  end

  #### show ####
  def show; puts "[#{@th.join(", ")}]"; end
end

########################################
# Mixture 分布
########################################
class Mixture
  #### new ####
  attr_reader :cp, :cc
  def initialize(_p, _c)
    @cp = Categorical.new(_p)
    @cc = Array.new(_c.size){|i| Categorical.new(_c[i]) }
  end

  #### sample ####
  def sample; @cc[ @cp.sample ].sample; end

  #### show ####
  def show
    @cp.show
    @cc.each do |c|
      c.show
    end
  end
end

########################################
# Markov Chain
########################################
class MarkovChain
  #### new ####
  attr_reader :cp, :cc, :s
  def initialize(_p, _c)
    @cp = Categorical.new(_p)
    @cc = Array.new(_c.size){|i| Categorical.new(_c[i]) }
    @s  = -1 # previous state
  end

  #### sample ####
  def sample; @s = @s < 0 ? @cp.sample : @cc[@s].sample; end

  #### reset ####
  def reset; @s = -1; end

  #### show ####
  def show
    puts "#{@s}"
    @cp.show
    @cc.each do |c|
      c.show
    end
  end
end
