#!/usr/bin/env ruby
# -*- coding: utf-8 -*-

########################################
# Gamma Distribution
########################################
class Gamma
  #### new ####
  def initialize(_k, _th)
    @k  = _k.to_f
    @th = _th.to_f
  end
  attr_reader :k, :th

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
  def show
    puts "#{@k}, #{@th}"
  end
end

########################################
# Gamma Distribution
########################################
class Dirichlet
  #### new ####
  def initialize(_al)
    @al = _al.clone
  end
  attr_reader :al

  #### sample ####
  def sample
    th = Array.new(@al.size){|i| Gamma.new(@al[i], 1.0).sample}
    sum = th.inject(:+)
    th.map!{|i| i / sum}
  end

  #### show ####
  def show
    puts "[#{@al.join(", ")}]"
  end
end

########################################
# Categorical 分布
########################################
class Categorical
  #### new ####
  def initialize(_theta)
    @th = _theta.clone

    #### 正規化 ####
    if (sum = @th.inject(:+)) != 1.0
      @th.map!{|i| i / sum}
    end
  end
  attr_reader :th

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
  def show
    puts "[#{@th.join(", ")}]"
  end
end

########################################
# Mixture 分布
########################################
class Mixture
  #### new ####
  def initialize(_p, _c)
    @cp = Categorical.new(_p)
    @cc = Array.new(_c.size){|i| Categorical.new(_c[i]) }
  end
  attr_reader :cp, :cc

  #### sample ####
  def sample
    k = @cp.sample
    @cc[k].sample
  end

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
  def initialize(_p, _c)
    @cp = Categorical.new(_p)
    @cc = Array.new(_c.size){|i| Categorical.new(_c[i]) }

    #### previous state ####
    @s = -1
  end
  attr_reader :cp, :cc, :s

  #### sample ####
  def sample
    if @s < 0
      t = @cp.sample
    else
      t = @cc[@s].sample
    end
    @s = t
    t
  end

  #### reset ####
  def reset
    @s = -1
  end

  #### show ####
  def show
    puts "#{@s}"
    @cp.show
    @cc.each do |c|
      c.show
    end
  end
end
