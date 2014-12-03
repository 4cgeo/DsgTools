<!DOCTYPE qgis PUBLIC 'http://mrcc.com/qgis.dtd' 'SYSTEM'><qgis version="2.6.0-Brighton" minimumScale="1" maximumScale="1" simplifyDrawingHints="0" minLabelScale="0" maxLabelScale="1e+08" simplifyDrawingTol="1" simplifyMaxScale="1" hasScaleBasedVisibilityFlag="0" simplifyLocal="1" scaleBasedLabelVisibilityFlag="0"> 
  <edittypes> 
     <edittype widgetv2type="TextEdit" name="OGC_FID"> 
      <widgetv2config IsMultiline="0" fieldEditable="0" UseHtml="0" labelOnTop="0"/> 
    </edittype> 
    <edittype widgetv2type="TextEdit" name="id"> 
      <widgetv2config IsMultiline="0" fieldEditable="0" UseHtml="0" labelOnTop="0"/> 
    </edittype>
    <edittype widgetv2type="TextEdit" name="id_area_verde"> 
      <widgetv2config IsMultiline="0" fieldEditable="0" UseHtml="0" labelOnTop="0"/> 
    </edittype>
    <edittype widgetv2type="ValueMap" name="geometriaaproximada">
      <widgetv2config fieldEditable="1" labelOnTop="0">
        <value key="N�o" value="0"/>
        <value key="Sim" value="1"/>
      </widgetv2config>
    </edittype> 
    <edittype widgetv2type="ValueMap" name="tipoveg">
      <widgetv2config fieldEditable="1" labelOnTop="0">
        <value key="Vegeta��o cultivada" value="1"/>
      </widgetv2config>
    </edittype> 
    <edittype widgetv2type="ValueMap" name="classificacaoporte">
      <widgetv2config fieldEditable="1" labelOnTop="0">
        <value key="Rasteira" value="2"/>
        <value key="Herb�cea" value="3"/>
        <value key="Arb�rea" value="4"/>
        <value key="Arbustiva" value="5"/>
        <value key="Desconhecida" value="95"/>
        <value key="Mista" value="97"/>
      </widgetv2config>
    </edittype> 
    <edittype widgetv2type="ValueMap" name="densidade">
      <widgetv2config fieldEditable="1" labelOnTop="0">
        <value key="Alta" value="1"/>
        <value key="Baixa" value="2"/>
        <value key="Desconhecida" value="95"/>
      </widgetv2config>
    </edittype> 
    <edittype widgetv2type="ValueMap" name="tipolavoura">
      <widgetv2config fieldEditable="1" labelOnTop="0">
        <value key="Anual" value="20"/>
        <value key="Perene" value="21"/>
        <value key="Semi-perene" value="22"/>
        <value key="Desconhecido" value="95"/>
      </widgetv2config>
    </edittype> 
    <edittype widgetv2type="ValueMap" name="terreno">
      <widgetv2config fieldEditable="1" labelOnTop="0">
        <value key="Seco" value="1"/>
        <value key="Irrigado" value="2"/>
        <value key="Inundado" value="3"/>
        <value key="Desconhecida" value="95"/>
      </widgetv2config>
    </edittype> 
    <edittype widgetv2type="ValueMap" name="finalidade">
      <widgetv2config fieldEditable="1" labelOnTop="0">
        <value key="Ornamental" value="2"/>
        <value key="Explora��o econ�mica" value="3"/>
        <value key="Subsist�ncia" value="4"/>
        <value key="Conserva��o ambiental" value="5"/>
        <value key="Desconhecida" value="95"/>
        <value key="Outros" value="99"/>
      </widgetv2config>
    </edittype> 
    <edittype widgetv2type="ValueMap" name="cultivopredominante">
      <widgetv2config fieldEditable="1" labelOnTop="0">
        <value key="Cultura rotativa" value="2"/>
        <value key="Milho" value="3"/>
        <value key="Banana" value="4"/>
        <value key="Laranja" value="5"/>
        <value key="Trigo" value="6"/>
        <value key="Abacate" value="7"/>
        <value key="Algod�o herb�ceo" value="8"/>
        <value key="Cana de a��car" value="9"/>
        <value key="Fumo" value="10"/>
        <value key="Soja" value="11"/>
        <value key="Batata inglesa" value="12"/>
        <value key="Mandioca" value="13"/>
        <value key="Feij�o" value="14"/>
        <value key="Arroz" value="15"/>
        <value key="Caf�" value="16"/>
        <value key="Cacau" value="17"/>
        <value key="Erva-mate" value="18"/>
        <value key="Palmeira" value="19"/>
        <value key="A�a�" value="20"/>
        <value key="Seringueira" value="21"/>
        <value key="Eucalipto" value="22"/>
        <value key="Ac�cia" value="23"/>
        <value key="Algaroba" value="24"/>
        <value key="Pinus" value="25"/>
        <value key="Pastagem cultivada" value="26"/>
        <value key="Hortali�as" value="27"/>
        <value key="Abacaxi ou anan�s" value="28"/>
        <value key="Arauc�ria" value="29"/>
        <value key="Carna�ba" value="30"/>
        <value key="Alfafa" value="31"/>
        <value key="Ma��" value="32"/>
        <value key="P�ssego" value="33"/>
        <value key="Juta" value="34"/>
        <value key="Cebola" value="35"/>
        <value key="Algod�o arb�reo" value="36"/>
        <value key="Alho" value="37"/>
        <value key="Amendoim" value="38"/>
        <value key="Aveia" value="39"/>
        <value key="Azeitona" value="40"/>
        <value key="Batata-doce" value="41"/>
        <value key="Caju" value="42"/>
        <value key="Centeio" value="43"/>
        <value key="Videira" value="44"/>
        <value key="Cevada" value="45"/>
        <value key="Ch�-da-�ndia" value="46"/>
        <value key="Coco-da-ba�a" value="47"/>
        <value key="Cravo da �ndia" value="48"/>
        <value key="Dend�" value="49"/>
        <value key="Ervilha" value="50"/>
        <value key="Fava" value="51"/>
        <value key="Figo" value="52"/>
        <value key="Flores" value="53"/>
        <value key="Girassol" value="54"/>
        <value key="Goiaba" value="55"/>
        <value key="Guaran�" value="56"/>
        <value key="Inhame" value="57"/>
        <value key="Lim�o" value="58"/>
        <value key="Linho" value="59"/>
        <value key="Malva" value="60"/>
        <value key="Mam�o" value="61"/>
        <value key="Mamona" value="62"/>
        <value key="Mandioca, aipim ou macaxeira" value="63"/>
        <value key="Manga" value="64"/>
        <value key="Maracuj�" value="65"/>
        <value key="Marmelo" value="66"/>
        <value key="Melancia" value="67"/>
        <value key="Mel�o" value="68"/>
        <value key="N�o identificado" value="96"/>
        <value key="Noz" value="70"/>
        <value key="Palmito" value="71"/>
        <value key="Pera" value="72"/>
        <value key="Pia�ava" value="73"/>
        <value key="Plantas ornamentais" value="74"/>
        <value key="Policultura" value="75"/>
        <value key="Rami" value="76"/>
        <value key="Sisal ou agave" value="77"/>
        <value key="Sorgo" value="78"/>
        <value key="Tangerina" value="79"/>
        <value key="Tomate" value="80"/>
        <value key="Triticale" value="81"/>
        <value key="Tungue" value="82"/>
        <value key="Urucum" value="83"/>
        <value key="Uva" value="84"/>
        <value key="Pimenta do reino" value="85"/>
        <value key="Outros" value="99"/>
      </widgetv2config>
    </edittype> 
  </edittypes>
</qgis>