
#
# class Fraction:
#      count = 0
#      def __init__(self,top,bottom):
#          self._num = top
#          self.den = bottom
#          Fraction.count +=1
#
#      @property
#      def num(self): #可通过 对象名.num访问成员变量_num，无法修改
#          return self._num
#
#      def gcd(self,m, n):
#          while m % n != 0:
#              oldm = m
#              oldn = n
#
#              m = oldn
#              n = oldm % oldn
#          return n
#
#      def __str__(self):
#          return str(self.num)+"/"+str(self.den)
#
#      def show(self):
#          print(self.num,"/",self.den)
#
#      def __add__(self,otherfraction):
#          newnum = self.num*otherfraction.den + \
#                       self.den*otherfraction.num
#          newden = self.den * otherfraction.den
#          common = self.gcd(newnum,newden)
#          return Fraction(newnum//common,newden//common)
#
#      def __eq__(self, other):
#          firstnum = self.num * other.den
#          secondnum = other.num * self.den
#
#          return firstnum == secondnum
# x = Fraction(1,2)
# y = Fraction(2,3)
# print(x+y)
# print(x == y)


class LogicGate:

    def __init__(self,n):
        self.name = n
        self.output = None

    def getName(self):
        return self.name

    def getOutput(self):
        self.output = self.performGateLogic()
        return self.output


class BinaryGate(LogicGate):

    def __init__(self,n):
        LogicGate.__init__(self,n)

        self.pinA = None
        self.pinB = None

    def getPinA(self):
        if self.pinA == None:
            return int(input("Enter Pin A input for gate "+self.getName()+"-->"))
        else:
            return self.pinA.getFrom().getOutput()

    def getPinB(self):
        if self.pinB == None:
            return int(input("Enter Pin B input for gate "+self.getName()+"-->"))
        else:
            return self.pinB.getFrom().getOutput()

    def setNextPin(self,source):
        if self.pinA == None:
            self.pinA = source
        else:
            if self.pinB == None:
                self.pinB = source
            else:
                print("Cannot Connect: NO EMPTY PINS on this gate")


class AndGate(BinaryGate):

    def __init__(self,n):
        BinaryGate.__init__(self,n)

    def performGateLogic(self):

        a = self.getPinA()
        b = self.getPinB()
        if a==1 and b==1:
            return 1
        else:
            return 0

class OrGate(BinaryGate):

    def __init__(self,n):
        BinaryGate.__init__(self,n)

    def performGateLogic(self):

        a = self.getPinA()
        b = self.getPinB()
        if a ==1 or b==1:
            return 1
        else:
            return 0

class UnaryGate(LogicGate):

    def __init__(self,n):
        LogicGate.__init__(self,n)

        self.pin = None

    def getPin(self):
        if self.pin == None:
            return int(input("Enter Pin input for gate "+self.getName()+"-->"))
        else:
            return self.pin.getFrom().getOutput()

    def setNextPin(self,source):
        if self.pin == None:
            self.pin = source
        else:
            print("Cannot Connect: NO EMPTY PINS on this gate")


class NotGate(UnaryGate):

    def __init__(self,n):
        UnaryGate.__init__(self,n)

    def performGateLogic(self):
        if self.getPin():
            return 0
        else:
            return 1


class Connector:

    def __init__(self, fgate, tgate):
        self.fromgate = fgate
        self.togate = tgate

        tgate.setNextPin(self)

    def getFrom(self):
        return self.fromgate

    def getTo(self):
        return self.togate


def main():
   g1 = AndGate("G1")
   g3 = OrGate("G3")
   c1 = Connector(g1,g3)
   print(c1.getTo().getOutput())

main()
