<!DOCTYPE qgis PUBLIC 'http://mrcc.com/qgis.dtd' 'SYSTEM'><qgis version="2.6.0-Brighton" minimumScale="1" maximumScale="1" simplifyDrawingHints="0" minLabelScale="0" maxLabelScale="1e+08" simplifyDrawingTol="1" simplifyMaxScale="1" hasScaleBasedVisibilityFlag="0" simplifyLocal="1" scaleBasedLabelVisibilityFlag="0"> 
  <edittypes> 
     <edittype widgetv2type="TextEdit" name="OGC_FID"> 
      <widgetv2config IsMultiline="0" fieldEditable="0" UseHtml="0" labelOnTop="0"/> 
    </edittype> 
    <edittype widgetv2type="TextEdit" name="id"> 
      <widgetv2config IsMultiline="0" fieldEditable="0" UseHtml="0" labelOnTop="0"/> 
    </edittype>
    <edittype widgetv2type="TextEdit" name="id_org_ext_mineral"> 
      <widgetv2config IsMultiline="0" fieldEditable="0" UseHtml="0" labelOnTop="0"/> 
    </edittype>
    <edittype widgetv2type="ValueMap" name="geometriaaproximada">
      <widgetv2config fieldEditable="1" labelOnTop="0">
        <value key="N�o" value="0"/>
        <value key="Sim" value="1"/>
      </widgetv2config>
    </edittype> 
    <edittype widgetv2type="ValueMap" name="tipoalterantrop">
      <widgetv2config fieldEditable="1" labelOnTop="0">
        <value key="Extrativismo mineral" value="27"/>
      </widgetv2config>
    </edittype> 
    <edittype widgetv2type="ValueMap" name="secaoativecon">
      <widgetv2config fieldEditable="1" labelOnTop="0">
        <value key="Ind�strias de transforma��o" value="1"/>
        <value key="Constru��o" value="2"/>
        <value key="Ind�strias extrativas" value="3"/>
        <value key="Desconhecido" value="95"/>
        <value key="Outros" value="99"/>
      </widgetv2config>
    </edittype> 
    <edittype widgetv2type="ValueMap" name="operacional">
      <widgetv2config fieldEditable="1" labelOnTop="0">
        <value key="N�o" value="0"/>
        <value key="Sim" value="1"/>
        <value key="Desconhecido" value="95"/>
      </widgetv2config>
    </edittype> 
    <edittype widgetv2type="ValueMap" name="situacaofisica">
      <widgetv2config fieldEditable="1" labelOnTop="0">
        <value key="Planejada" value="1"/>
        <value key="Constru�da" value="2"/>
        <value key="Abandonada" value="3"/>
        <value key="Destru�da" value="4"/>
        <value key="Em constru��o" value="5"/>
        <value key="Constru�da, mas em obras" value="6"/>
        <value key="Desconhecida" value="95"/>
        <value key="N�o aplic�vel" value="97"/>
      </widgetv2config>
    </edittype> 
    <edittype widgetv2type="ValueMap" name="tipoextmin">
      <widgetv2config fieldEditable="1" labelOnTop="0">
        <value key="Po�o de petr�leo" value="1"/>
        <value key="Mina" value="2"/>
        <value key="Po�o para �gua subterr�nea" value="3"/>
        <value key="Salina" value="4"/>
        <value key="Garimpo" value="5"/>
        <value key="Desconhecido" value="95"/>
        <value key="Outros" value="99"/>
      </widgetv2config>
    </edittype> 
    <edittype widgetv2type="ValueMap" name="tipoproduto">
      <widgetv2config fieldEditable="1" labelOnTop="0">
        <value key="Sal-gema" value="1"/>
        <value key="Terras raras" value="2"/>
        <value key="Tit�nio" value="3"/>
        <value key="Top�zio" value="4"/>
        <value key="Tungst�nio" value="5"/>
        <value key="Turmalina" value="6"/>
        <value key="T�rio" value="7"/>
        <value key="Ur�nio" value="8"/>
        <value key="Opala" value="9"/>
        <value key="Zinco" value="10"/>
        <value key="Zirc�nio" value="11"/>
        <value key="N�quel" value="12"/>
        <value key="Querosene" value="13"/>
        <value key="�gua mineral" value="14"/>
        <value key="�leo diesel" value="15"/>
        <value key="Vermiculita" value="16"/>
        <value key="�gata" value="17"/>
        <value key="�gua" value="18"/>
        <value key="Ni�bio" value="19"/>
        <value key="Rocha ornamental" value="20"/>
        <value key="Ouro" value="24"/>
        <value key="Petr�leo" value="25"/>
        <value key="Pedra preciosa" value="26"/>
        <value key="G�s" value="27"/>
        <value key="Gr�o" value="28"/>
        <value key="Alexandrita" value="29"/>
        <value key="Ametista" value="30"/>
        <value key="Amianto" value="31"/>
        <value key="Argila" value="32"/>
        <value key="Barita" value="33"/>
        <value key="Bentonita" value="34"/>
        <value key="Calc�rio" value="35"/>
        <value key="Carv�o vegetal" value="36"/>
        <value key="Caulim" value="37"/>
        <value key="Vinhoto" value="38"/>
        <value key="Estrume" value="39"/>
        <value key="Cascalho" value="40"/>
        <value key="Chumbo" value="41"/>
        <value key="Inseticida" value="42"/>
        <value key="Folhagem" value="43"/>
        <value key="�gua marinha" value="44"/>
        <value key="Pedra (brita)" value="45"/>
        <value key="Granito" value="46"/>
        <value key="M�rmore" value="47"/>
        <value key="Bauxita" value="48"/>
        <value key="Mangan�s" value="49"/>
        <value key="Talco" value="50"/>
        <value key="Chorume" value="51"/>
        <value key="Gasolina" value="52"/>
        <value key="�lcool" value="53"/>
        <value key="Citrino" value="54"/>
        <value key="Cobre" value="55"/>
        <value key="Carv�o mineral" value="56"/>
        <value key="Sal" value="57"/>
        <value key="Turfa" value="58"/>
        <value key="Esc�ria" value="59"/>
        <value key="Ferro" value="60"/>
        <value key="Crisoberilo" value="61"/>
        <value key="Prata" value="62"/>
        <value key="Cristal de rocha" value="63"/>
        <value key="Forragem" value="64"/>
        <value key="Saibro/pi�arra" value="65"/>
        <value key="Areia" value="66"/>
        <value key="Cromo" value="67"/>
        <value key="Diamante" value="68"/>
        <value key="Diatomita" value="69"/>
        <value key="Dolomito" value="70"/>
        <value key="Esgoto" value="71"/>
        <value key="Esmeralda" value="72"/>
        <value key="Estanho" value="73"/>
        <value key="Feldspato" value="74"/>
        <value key="Fosfato" value="75"/>
        <value key="Gipsita" value="76"/>
        <value key="Grafita" value="77"/>
        <value key="Granada" value="78"/>
        <value key="Lixo domiciliar e comercial" value="79"/>
        <value key="Lixo s�ptico" value="80"/>
        <value key="Lixo t�xico" value="81"/>
        <value key="L�tio" value="82"/>
        <value key="Magnesita" value="83"/>
        <value key="Mica" value="84"/>
        <value key="Desconhecido" value="95"/>
        <value key="Outros" value="99"/>
      </widgetv2config>
    </edittype> 
    <edittype widgetv2type="ValueMap" name="tipopocomina">
      <widgetv2config fieldEditable="1" labelOnTop="0">
        <value key="Horizontal" value="3"/>
        <value key="Vertical" value="4"/>
        <value key="Desconhecido" value="95"/>
        <value key="N�o aplic�vel" value="97"/>
      </widgetv2config>
    </edittype> 
    <edittype widgetv2type="ValueMap" name="procextracao">
      <widgetv2config fieldEditable="1" labelOnTop="0">
        <value key="Manual" value="1"/>
        <value key="Mecanizado" value="2"/>
        <value key="Desconhecido" value="95"/>
      </widgetv2config>
    </edittype> 
    <edittype widgetv2type="ValueMap" name="atividade">
      <widgetv2config fieldEditable="1" labelOnTop="0">
        <value key="Prospec��o" value="9"/>
        <value key="Produ��o" value="10"/>
        <value key="Desconhecida" value="95"/>
      </widgetv2config>
    </edittype> 
    <edittype widgetv2type="ValueMap" name="formaextracao">
      <widgetv2config fieldEditable="1" labelOnTop="0">
        <value key="A c�u aberto" value="1"/>
        <value key="Subterr�nea" value="2"/>
        <value key="Desconhecida" value="95"/>
      </widgetv2config>
    </edittype> 
  </edittypes>
</qgis>