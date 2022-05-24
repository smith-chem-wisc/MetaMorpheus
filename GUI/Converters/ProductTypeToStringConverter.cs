using Newtonsoft.Json;
using Newtonsoft.Json.Linq;
using OxyPlot;
using Proteomics.Fragmentation;
using System;
using System.Collections;
using System.Collections.Generic;
using System.Diagnostics;
using System.Globalization;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace MetaMorpheusGUI
{
    public class ProductTypeToStringConverter : BaseValueConverter<ProductTypeToStringConverter>
    {
        public override object Convert(object value, Type targetType, object parameter, CultureInfo culture)
        {
            switch (((KeyValuePair<ProductType, OxyColor>)value).Key) 
            {
                case ProductType.a:
                    return new string("a");
                case ProductType.aStar:
                    return new string("a Star");
                case ProductType.aDegree:
                    return new string("a Degree");
                case ProductType.b:
                    return new string("b");
                case ProductType.bStar:
                    return new string("b Star");
                case ProductType.bDegree:
                    return new string("b Degree");
                case ProductType.c:
                    return new string("c");
                case ProductType.x:
                    return new string("x");
                case ProductType.y:
                    return new string("y");
                case ProductType.yStar:
                    return new string("y Star");
                case ProductType.yDegree:
                    return new string("y Degree");
                case ProductType.zPlusOne:
                    return new string("z Plus One");
                case ProductType.zDot:
                    return new string("z Dot");
                case ProductType.M:
                    return new string("M");
                case ProductType.D:
                    return new string("D");
                case ProductType.Ycore:
                    return new string("Y core");
                case ProductType.Y:
                    return new string("Y");
                default:
                    Debugger.Break();
                    return null;
            }
        }

        public override object ConvertBack(object value, Type targetType, object parameter, CultureInfo culture)
        {
            throw new NotImplementedException();
        }
    }
}
